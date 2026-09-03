#pragma once

#include <vector>
#include <utility>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "diskpp/methods/hho"
#include "diskpp/solvers/direct_solvers.hpp"

#include "spatial_operator.hpp"
#include "contact_composer.hpp"
#include "hho_params.hpp"

namespace hho_contact {

template<typename T>
struct hho_state
{
    std::vector<disk::dynamic_vector<T>> cells;
    disk::dynamic_vector<T>              faces;
};

/* Per-solve iteration record, filled by solve().  For a Newton-only
 * problem picard_iters = 1 and newton_per_picard has a single entry.     */
struct solve_stats
{
    bool                converged   = false;
    size_t              picard_iters = 0;   // outer fixed-point passes
    size_t              newton_total = 0;   // sum over all passes
    std::vector<size_t> newton_per_picard;  // one entry per Picard pass
};

template<typename Mesh>
class NewtonSolver
{
public:
    using T           = typename Mesh::coordinate_type;
    using matrix_type = disk::dynamic_matrix<T>;
    using vector_type = disk::dynamic_vector<T>;
    using cell_type   = typename Mesh::cell_type;
    using bc_type     = disk::vector_boundary_conditions<Mesh>;

    using assembler_type =
        decltype(disk::make_vector_primal_hho_assembler(
                     std::declval<const Mesh&>(),
                     std::declval<const disk::MeshDegreeInfo<Mesh>&>(),
                     std::declval<const bc_type&>()));

private:
    const spatial_operator<Mesh>& op;
    const ContactTerms<Mesh>&     ct;
    const bc_type&                bnd;        // state gather (loads, faces)
    const bc_type&                bnd_incr;   // increment assembly / scatter

    assembler_type assembler;

    std::vector<matrix_type> A;    // K - consistency + linear penalties   (cached)
    std::vector<vector_type> L;    // body + Neumann loads                 (cached)
    std::vector<matrix_type> AL;   // per-iteration condensation storage
    std::vector<vector_type> bL;

public:
    template<typename BodyF>
    NewtonSolver(const spatial_operator<Mesh>& op_,
                 const ContactTerms<Mesh>&     ct_,
                 const bc_type& bnd_state, const bc_type& bnd_increment,
                 BodyF&& body)
        : op(op_), ct(ct_), bnd(bnd_state), bnd_incr(bnd_increment),
          assembler(disk::make_vector_primal_hho_assembler(op_.msh, op_.di,
                                                           bnd_state))
    {
        const size_t nc = op.num_cells();
        A.resize(nc); L.resize(nc); AL.resize(nc); bL.resize(nc);

        size_t ci = 0;
        for (auto& cl : op.msh)
        {
            // ---- constant system matrix: every linear term leaves the loop
            A[ci] = op.K[ci]
                    - ct.consistency(cl, ci)
                    + ct.linear_jacobian(cl, ci);

            // ---- loads
            auto cb = disk::make_vector_monomial_basis(op.msh, cl, op.cell_deg);
            L[ci]   = vector_type::Zero(op.total_dofs(cl));
            L[ci].head(op.cbs) = disk::make_rhs(op.msh, cl, cb, body);

            size_t ofs = op.cbs;
            for (auto& fc : disk::faces(op.msh, cl))
            {
                const auto   fdi = op.di.degreeInfo(op.msh, fc);
                const size_t nfd = disk::vector_face_dofs(op.msh, fdi);
                if (bnd.is_neumann_face(fc))
                {
                    auto fb = disk::make_vector_monomial_basis(op.msh, fc,
                                                               fdi.degree());
                    L[ci].segment(ofs, nfd) += disk::make_rhs(op.msh, fc, fb,
                        bnd.neumann_boundary_func(op.msh.lookup(fc)));
                }
                ofs += nfd;
            }
            ci++;
        }
    }

    // ---- state ----------------------------------------------------------
    hho_state<T> make_state() const
    {
        hho_state<T> st;
        st.cells.assign(op.num_cells(), vector_type::Zero(op.cbs));
        st.faces = vector_type::Zero(assembler.LHS.rows());
        return st;
    }

    vector_type gather(const cell_type& cl, size_t ci,
                       const hho_state<T>& st) const
    {
        vector_type u = vector_type::Zero(A[ci].rows());
        u.head(op.cbs) = st.cells[ci];
        vector_type uF = assembler.take_local_solution(op.msh, cl, bnd, st.faces);
        u.tail(uF.size()) = uF;
        return u;
    }

    // ---- driver entry point ---------------------------------------------
    bool solve(hho_state<T>& st, const solver_opts& o)
    {
        solve_stats stats;
        return solve(st, o, stats);
    }

    bool solve(hho_state<T>& st, const solver_opts& o, solve_stats& stats)
    {
        stats = solve_stats{};

        if (!ct.needs_picard())
        {
            size_t nit = 0;
            const bool ok = newton_loop(st, nullptr, o, nit, o.verbose);
            stats.converged    = ok;
            stats.picard_iters = 1;
            stats.newton_total = nit;
            stats.newton_per_picard = { nit };
            if (o.verbose)
            {
                if (!ok) std::cout << "  *** NEWTON DID NOT CONVERGE ***\n";
                std::cout << "  Newton iterations: " << nit << "\n";
            }
            return ok;
        }

        // ---- Picard on the frozen Coulomb threshold ----------------------
        hho_state<T> frozen = st;
        size_t total_newton = 0;
        bool   picard_ok    = false;

        for (size_t pit = 0; pit < o.max_picard; pit++)
        {
            size_t nit = 0;
            const bool newton_ok = newton_loop(st, &frozen, o, nit, false);
            total_newton += nit;
            stats.newton_per_picard.push_back(nit);

            T diff2 = (st.faces - frozen.faces).squaredNorm();
            for (size_t ci = 0; ci < st.cells.size(); ci++)
                diff2 += (st.cells[ci] - frozen.cells[ci]).squaredNorm();
            const T pdiff = std::sqrt(diff2);

            if (o.verbose)
                std::cout << "  picard " << std::setw(2) << pit
                          << "   newton " << std::setw(2) << nit
                          << (newton_ok ? " " : "*")
                          << "   |u - u_picard| = " << std::scientific
                          << std::setprecision(3) << pdiff << "\n";

            frozen = st;

            if (pdiff < o.tol_picard) { picard_ok = true; break; }
        }

        stats.converged    = picard_ok;
        stats.picard_iters = stats.newton_per_picard.size();
        stats.newton_total = total_newton;

        if (o.verbose)
        {
            if (!picard_ok) std::cout << "  *** PICARD DID NOT CONVERGE ***\n";
            std::cout << "  total inner Newton iterations: "
                      << total_newton << "\n";
        }
        return picard_ok;
    }

private:
    /* One condensed Newton step; returns (|du_faces|, |u_faces|).
     * Rhist is the sliding window of accepted-iterate residual norms used
     * by the non-monotone line search below; it lives in the caller
     * (newton_loop), reset at the start of each Newton loop.             */
    std::pair<T,T> step(hho_state<T>& st, const hho_state<T>* frozen,
                        std::vector<T>& Rhist)
    {
        const bool nonlin = ct.has_nonlinear();

        // ---- pass 1: uncondensed local tangent Jbase[ci] and residual
        // R_loc[ci] at the CURRENT state. These don't change while we try
        // different damping/step-size combinations below -- only the
        // condensation (which sees Jbase + mu*I) and the global solve do.
        std::vector<matrix_type> Jbase(op.num_cells());
        std::vector<vector_type> Rloc(op.num_cells());
        T R0sq = 0;

        size_t ci = 0;
        for (auto& cl : op.msh)
        {
            vector_type u_loc = gather(cl, ci, st);

            matrix_type J = A[ci];
            vector_type R = A[ci]*u_loc - L[ci];

            if (nonlin)
            {
                // Coulomb reads p; everyone else ignores it, so avoid a
                // second gather unless an outer iteration is running.
                if (frozen)
                {
                    vector_type p_loc = gather(cl, ci, *frozen);
                    ct.accumulate_nonlinear(cl, ci, u_loc, p_loc, J, R);

// --- FD check: ¿J es la derivada de R? (1a celda de contacto, 1a vez) ---
{
    static bool fd_done = false;
    if (!fd_done && (J - A[ci]).norm() > 1e-12)   // celda con termino de contacto
    {
        fd_done = true;
        const double eps = 1e-6;
        vector_type d = vector_type::Random(u_loc.size());
        d /= d.norm();
        vector_type u2 = u_loc + eps*d;            // materializar
        matrix_type J2 = A[ci];
        vector_type R2 = A[ci]*u2 - L[ci];
        ct.accumulate_nonlinear(cl, ci, u2, p_loc, J2, R2);   // orden correcto
        vector_type fd = (R2 - R)/eps;
        vector_type an = J*d;
        std::cout << "  [dbg-fd] |fd - J*d|/|J*d| = "
                  << (fd-an).norm()/std::max(an.norm(),1e-30) << "\n";
    }
}


                }
                else
                    ct.accumulate_nonlinear(cl, ci, u_loc, u_loc, J, R);
            }

            // --- diagnostico temporal: definitud local, solo primera pasada ---
            {
                static int probe_budget = 280;   // = nº de celdas -> solo 1ª iteración
                if (probe_budget > 0) {
                    probe_budget--;
                    Eigen::SelfAdjointEigenSolver<matrix_type> es(0.5*(J + J.transpose()));
                    const double mn = es.eigenvalues().minCoeff();
                    const double mx = es.eigenvalues().maxCoeff();
                    if (mn < -1e-10 * std::abs(mx))
                        std::cout << "  [dbg-eig] ci=" << ci << "  min=" << mn
                                << "  max=" << mx << "\n";
                }
            }

            R0sq += R.squaredNorm();
            Jbase[ci] = J;
            Rloc[ci]  = R;
            ci++;
        }

        const T R0 = std::sqrt(R0sq);
        const size_t window = 5;
        const T Rmax = Rhist.empty()
                     ? R0 : *std::max_element(Rhist.begin(), Rhist.end());
        Rhist.push_back(R0);
        if (Rhist.size() > window) Rhist.erase(Rhist.begin());

        // ---- pass 2: damped semismooth Newton.
        //
        // A plain backtracking line search can only rescale the Newton
        // direction d = J^-1(-R) by a scalar alpha; it cannot fix a d that
        // simply is not a descent direction for the residual at all, which
        // is what happens once the contact active set (cell variant) has
        // settled somewhere the linearised model stops predicting the true
        // nonlinear residual well. Levenberg-Marquardt-style damping --
        // solving (Jbase + mu*diag_scale*I) d = -R instead -- rotates the
        // direction towards steepest descent as mu grows, which for large
        // enough mu is *guaranteed* to be a descent direction, unlike any
        // rescaling of the undamped direction. mu escalates only when the
        // (cheap) backtracking line search fails outright at the current mu.
        hho_state<T> best        = st;
        T            best_R      = R0;
        bool         accepted    = false;
        T            accepted_alpha = T(1);
        T            accepted_mu    = T(0);

        T mu = T(0);
        const size_t max_mu_tries = 6;

        for (size_t mu_try = 0; mu_try <= max_mu_tries && !accepted; mu_try++)
        {
            assembler.initialize();

            ci = 0;
            for (auto& cl : op.msh)
            {
                matrix_type Jreg = Jbase[ci];
                if (mu > 0)
                {
                    const T dscale = std::max(Jbase[ci].diagonal().cwiseAbs().maxCoeff(),
                                              T(1e-30));
                    Jreg += (mu*dscale) * matrix_type::Identity(Jreg.rows(), Jreg.cols());
                }

                vector_type negR = -Rloc[ci];
                auto scT = disk::make_vector_static_condensation_withMatrix(
                               op.msh, cl, op.di, Jreg, negR);
                auto sc  = std::get<0>(scT);
                AL[ci]   = std::get<1>(scT);
                bL[ci]   = std::get<2>(scT);

                assembler.assemble(op.msh, cl, bnd_incr, sc.first, sc.second);
                ci++;
            }
            assembler.finalize();

            vector_type du_faces = vector_type::Zero(assembler.LHS.rows());
#ifdef HAVE_PARDISO
            disk::solvers::sparse_lu(assembler.LHS, assembler.RHS, du_faces,
                                     disk::solvers::direct_solver::pardiso);
#else
            disk::solvers::sparse_lu(assembler.LHS, assembler.RHS, du_faces,
                                     disk::solvers::direct_solver::sparselu);
#endif

            std::vector<vector_type> dcells(op.num_cells());
            ci = 0;
            for (auto& cl : op.msh)
            {
                vector_type duF = assembler.take_local_solution(op.msh, cl,
                                                                 bnd_incr, du_faces);
                dcells[ci] = bL[ci] - AL[ci]*duF;
                ci++;
            }

            // cheap backtracking on this (fixed-mu) direction
            const size_t max_ls = 5;
            T alpha = T(1);
            for (size_t ls = 0; ls <= max_ls; ls++)
            {
                hho_state<T> trial = st;
                trial.faces += alpha * du_faces;
                for (size_t c = 0; c < op.num_cells(); c++)
                    trial.cells[c] += alpha * dcells[c];

                const T Rtrial = residual_norm(trial, frozen);

                if (Rtrial < best_R) { best = trial; best_R = Rtrial; }

                if (Rtrial <= Rmax)
                {
                    accepted       = true;
                    accepted_alpha = alpha;
                    accepted_mu    = mu;
                    break;
                }
                if (ls == max_ls) break;
                alpha *= T(0.5);
            }

            mu = (mu == T(0)) ? T(1e-4) : mu*T(10);
        }

        // Fall back to the best trial found across every (mu, alpha) tried,
        // even if it never satisfied the non-monotone criterion: guarantees
        // the state keeps moving instead of freezing at st.
        const hho_state<T> old_faces_state = st;
        st = best;

        if (!accepted || accepted_mu > 0 || accepted_alpha < T(1))
            std::cout << "    damped newton: mu = " << accepted_mu
                      << "   alpha = " << accepted_alpha
                      << "   accepted = " << accepted
                      << "   |R0| = " << R0 << "   |Rmax| = " << Rmax
                      << "   |R| = " << best_R << "\n";

        vector_type dfaces = st.faces - old_faces_state.faces;
        return { dfaces.norm(), st.faces.norm() };
    }

    /* Global residual norm (uncondensed, summed over cells) at a trial
     * state -- the merit function the line search backtracks on.        */
    T residual_norm(const hho_state<T>& st, const hho_state<T>* frozen) const
    {
        T s = 0;
        size_t ci = 0;
        for (auto& cl : op.msh)
        {
            vector_type u_loc = gather(cl, ci, st);

            matrix_type J = A[ci];
            vector_type R = A[ci]*u_loc - L[ci];

            if (ct.has_nonlinear())
            {
                if (frozen)
                {
                    vector_type p_loc = gather(cl, ci, *frozen);
                    ct.accumulate_nonlinear(cl, ci, u_loc, p_loc, J, R);
                }
                else
                    ct.accumulate_nonlinear(cl, ci, u_loc, u_loc, J, R);
            }

            s += R.squaredNorm();
            ci++;
        }
        return std::sqrt(s);
    }

    bool newton_loop(hho_state<T>& st, const hho_state<T>* frozen,
                     const solver_opts& o, size_t& nit, bool trace)
    {
        std::vector<T> Rhist;   // non-monotone line search window, fresh per loop

        for (nit = 0; nit < o.max_newton; nit++)
        {
            auto [du, u] = step(st, frozen, Rhist);

            if (trace)
                std::cout << "  newton " << std::setw(2) << nit
                          << "   |du| = " << std::scientific
                          << std::setprecision(3) << du << "\n";

            // absolute OR relative: whichever is looser at the current scale
            if (du < o.tol_newton || du < o.tol_newton_rel*(u + 1e-14))
            { nit++; return true; }
        }
        return false;
    }
};

} // namespace hho_contact