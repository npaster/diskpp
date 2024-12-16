/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

// NewtonRaphson iteration

#pragma once

#include <iostream>
#include <math.h>
#include <string>
#include <vector>

#include "diskpp/adaptivity/adaptivity.hpp"
#include "diskpp/bases/bases.hpp"
#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolverComput.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolverDynamic.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolverInformations.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolverParameters.hpp"
#include "diskpp/mechanics/NewtonSolver/StabilizationManager.hpp"
#include "diskpp/mechanics/NewtonSolver/TimeManager.hpp"
#include "diskpp/mechanics/NewtonSolver/Fields.hpp"
#include "diskpp/mechanics/behaviors/laws/behaviorlaws.hpp"
#include "diskpp/methods/hho"
#include "diskpp/solvers/solver.hpp"

#include "mumps.hpp"


namespace disk
{

namespace mechanics
{

/**
 * @brief Newton-Raphson iteration for nonlinear solid mechanics
 *
 *  Specialized for HHO methods
 *
 *  Options :  - small and finite deformations
 *             - plasticity, hyperelasticity (various laws)
 *
 * @tparam MeshType type of the mesh
 */
template <typename MeshType>
class NewtonIteration
{
    typedef MeshType mesh_type;
    typedef typename mesh_type::cell cell_type;
    typedef typename mesh_type::coordinate_type scalar_type;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    typedef NewtonSolverParameter<scalar_type> param_type;
    typedef vector_boundary_conditions<mesh_type> bnd_type;
    typedef Behavior<mesh_type> behavior_type;

    typedef vector_mechanics_hho_assembler<mesh_type> assembler_type;
    typedef mechanical_computation<mesh_type> elem_type;
    typedef dynamic_computation<mesh_type> dyna_type;

    vector_type m_system_displ;

    assembler_type m_assembler;

    std::vector<vector_type> m_bL;
    std::vector<matrix_type> m_AL;

    TimeStep<scalar_type> m_time_step;

    dyna_type m_dyna;

    scalar_type m_F_int, m_sol_norm;

    bool m_verbose;

    matrix_type
    _gradrec(const mesh_type &msh,
             const cell_type &cl,
             const param_type &rp,
             const MeshDegreeInfo<mesh_type> &degree_infos,
             const bool &small_def,
             const std::vector<matrix_type> &gradient_precomputed) const
    {
        if (rp.m_precomputation)
        {
            const auto cell_i = msh.lookup(cl);
            return gradient_precomputed[cell_i];
        }

        if (small_def)
        {
            const auto gradrec_sym_full = make_matrix_hho_symmetric_gradrec(msh, cl, degree_infos);
            return gradrec_sym_full.first;
        }
        else
        {
            const auto gradrec_full = make_matrix_hho_gradrec(msh, cl, degree_infos);
            return gradrec_full.first;
        }

        return matrix_type();
    }

    matrix_type
    _stab(const mesh_type &msh,
          const cell_type &cl,
          const param_type &rp,
          const MeshDegreeInfo<mesh_type> &degree_infos,
          const std::vector<matrix_type> &stab_precomputed) const
    {
        if (rp.m_precomputation)
        {
            const auto cell_i = msh.lookup(cl);
            return stab_precomputed.at(cell_i);
        }
        else
        {
            switch (rp.m_stab_type)
            {
            case StabilizationType::HHO_SYM: {
                const auto recons = make_vector_hho_symmetric_laplacian(msh, cl, degree_infos);
                return make_vector_hho_stabilization(msh, cl, recons.first,
                                                     degree_infos);
                break;
            }
            case StabilizationType::HHO: {
                const auto recons_scalar = make_scalar_hho_laplacian(msh, cl, degree_infos);
                return make_vector_hho_stabilization_optim(msh, cl, recons_scalar.first, degree_infos);
                break;
            }
            case StabilizationType::HDG: {
                return make_vector_hdg_stabilization(msh, cl, degree_infos);
                break;
            }
            case StabilizationType::DG: {
                return make_vector_dg_stabilization(msh, cl, degree_infos);
                break;
            }
            case StabilizationType::NO: {
                break;
            }
            default:
                throw std::invalid_argument("Unknown stabilization");
            }
        }

        return matrix_type();
    }

public:
    NewtonIteration(const mesh_type &msh,
                    const bnd_type &bnd,
                    const param_type &rp,
                    const MeshDegreeInfo<mesh_type> &degree_infos,
                    const TimeStep<scalar_type> &current_step) : m_verbose(rp.m_verbose),
                                                                 m_time_step(current_step),
                                                                 m_dyna(rp)
    {
        m_AL.clear();
        m_AL.resize(msh.cells_size());

        m_bL.clear();
        m_bL.resize(msh.cells_size());

        m_assembler = assembler_type(msh, degree_infos, bnd);
    }

    bool
    verbose(void) const
    {
        return m_verbose;
    }

    void
    verbose(bool v)
    {
        m_verbose = v;
    }

    void initialize(const mesh_type &msh, MultiTimeField<scalar_type> &fields) {
        m_dyna.prediction(msh, m_time_step, fields);
    }

    template <typename LoadFunction>
    AssemblyInfo compute_cells(const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                               const MeshDegreeInfo<mesh_type> &degree_infos,
                               const LoadFunction &lf,
                               const std::vector<matrix_type> &gradient_precomputed,
                               const std::vector<matrix_type> &stab_precomputed,
                               behavior_type &behavior, StabCoeffManager<scalar_type> &stab_manager,
                               MultiTimeField<scalar_type> &fields) {
        elem_type elem;
        AssemblyInfo ai;

        if (m_dyna.isExplicit()) {

            const bool small_def = (behavior.getDeformation() == SMALL_DEF);

            // Like if it is an implicit scheme
            const auto previous_time = m_time_step.start_time();
            const auto dt = m_time_step.increment_time();

            const auto depl_prev = fields.getField(-1, FieldName::DEPL);
            const auto deplT_prev = fields.getField(-1, FieldName::DEPL_CELLS);
            const auto vite_prev = fields.getField(-1, FieldName::VITE_CELLS);

            auto depl_curr = fields.getCurrentField(FieldName::DEPL);

            std::vector<vector_type> vite, acce, depl;
            vite.reserve(msh.cells_size());
            acce.reserve(msh.cells_size());
            depl.reserve(msh.cells_size());

            auto rlf = [&lf,
                        previous_time](const point<scalar_type, mesh_type::dimension> &p) -> auto {
                return lf(p, previous_time);
            };

            timecounter tc, ttot;

            ttot.tic();

            for (auto &cl : msh) {
                const auto cell_i = msh.lookup(cl);

                const auto huT = depl_prev.at(cell_i);

                // Gradient Reconstruction
                // std::cout << "Grad" << std::endl;
                tc.tic();
                matrix_type GT =
                    _gradrec(msh, cl, rp, degree_infos, small_def, gradient_precomputed);
                tc.toc();
                ai.m_time_gradrec += tc.elapsed();

                // Mechanical Computation

                tc.tic();
                // std::cout << "Elem" << std::endl;
                elem.compute(msh, cl, bnd, rp, degree_infos, rlf, GT, huT, m_time_step, behavior,
                             stab_manager, small_def, false);

                vector_type rhs_all = elem.RTF;

                tc.toc();
                ai.m_time_elem += tc.elapsed();
                ai.m_time_law += elem.time_law;

                // Stabilisation Contribution
                // std::cout << "Stab" << std::endl;
                tc.tic();
                if (rp.m_stab) {
                    const auto beta_s = stab_manager.getValue(msh, cl);

                    matrix_type stab = beta_s * _stab(msh, cl, rp, degree_infos, stab_precomputed);

                    // std::cout << beta_s << std::endl;

                    rhs_all -= stab * huT;
                }
                tc.toc();
                ai.m_time_stab += tc.elapsed();

                matrix_type lhs = (1.0 / (dt * dt)) * m_dyna.mass_matrix(msh, cl, degree_infos);

                const auto num_cell_dofs = lhs.rows();

                const vector_type depl_pred = deplT_prev.at(cell_i) + dt * vite_prev.at(cell_i);

                const vector_type rhs = rhs_all.head(num_cell_dofs) + lhs * depl_pred;

                const vector_type uT = lhs.ldlt().solve(rhs);
                const vector_type aT = (1.0 / (dt * dt)) * (uT - depl_pred);
                const vector_type vT = vite_prev.at(cell_i) + dt * aT;

                acce.push_back(aT);
                vite.push_back(vT);
                depl.push_back(uT);

                depl_curr[cell_i].head(num_cell_dofs) = uT;
            }

            fields.setCurrentField(FieldName::VITE_CELLS, vite);
            fields.setCurrentField(FieldName::ACCE_CELLS, acce);
            fields.setCurrentField(FieldName::DEPL_CELLS, depl);
            fields.setCurrentField(FieldName::DEPL, depl_curr);

            ttot.toc();
            ai.m_time_assembly = ttot.elapsed();
        }
        return ai;
    }

    template <typename LoadFunction>
    AssemblyInfo
    assemble(const mesh_type &msh,
             const bnd_type &bnd,
             const param_type &rp,
             const MeshDegreeInfo<mesh_type> &degree_infos,
             const LoadFunction &lf,
             const std::vector<matrix_type> &gradient_precomputed,
             const std::vector<matrix_type> &stab_precomputed,
             behavior_type &behavior,
             StabCoeffManager<scalar_type> &stab_manager,
             const MultiTimeField<scalar_type> &fields)
    {
        elem_type elem;
        AssemblyInfo ai;

        // set RHS to zero
        m_assembler.initialize();
        m_F_int = 0.0;

        const bool small_def = (behavior.getDeformation() == SMALL_DEF);

        // Like if it is an implicit scheme
        auto current_time = m_time_step.end_time();
        auto depl = fields.getCurrentField(FieldName::DEPL);
        auto depl_faces = fields.getCurrentField(FieldName::DEPL_FACES);

        std::vector<vector_type> acce_cells;
        if (m_dyna.enable()) {
            acce_cells = fields.getCurrentField(FieldName::ACCE_CELLS);
        }

        m_sol_norm = 0.0;
        for (auto &uF : depl_faces)
        {
            m_sol_norm += uF.squaredNorm();
        }

        auto rlf = [&lf, &current_time](const point<scalar_type, mesh_type::dimension> &p) -> auto {
            return lf(p, current_time);
        };

        timecounter tc, ttot;

        ttot.tic();

        for (auto &cl : msh)
        {
            const auto cell_i = msh.lookup(cl);

            const auto huT = depl.at(cell_i);

            // Gradient Reconstruction
            // std::cout << "Grad" << std::endl;
            tc.tic();
            matrix_type GT = _gradrec(msh, cl, rp, degree_infos, small_def, gradient_precomputed);
            tc.toc();
            ai.m_time_gradrec += tc.elapsed();

            // Mechanical Computation

            tc.tic();
            // std::cout << "Elem" << std::endl;
            elem.compute(msh, cl, bnd, rp, degree_infos, rlf, GT, huT, m_time_step, behavior,
                         stab_manager, small_def);

            matrix_type lhs = elem.K_int;
            vector_type rhs = elem.RTF;
            m_F_int += elem.F_int.squaredNorm();

            tc.toc();
            ai.m_time_elem += tc.elapsed();
            ai.m_time_law += elem.time_law;

            // Stabilisation Contribution
            // std::cout << "Stab" << std::endl;
            tc.tic();
            if (rp.m_stab)
            {
                const auto beta_s = stab_manager.getValue(msh, cl);

                matrix_type stab = beta_s * _stab(msh, cl, rp, degree_infos, stab_precomputed);

                // std::cout << beta_s << std::endl;

                lhs += stab;
                rhs -= stab * huT;
            }
            tc.toc();
            ai.m_time_stab += tc.elapsed();

            // Dynamic contribution
            if (m_dyna.enable()) {
                if (!m_dyna.isExplicit()) {
                    m_dyna.compute(msh, cl, degree_infos, huT, acce_cells.at(cell_i), m_time_step);
                    lhs += m_dyna.K_iner;
                    rhs += m_dyna.R_iner;
                    ai.m_time_dyna += m_dyna.time_dyna;
                } else {
                    const auto num_cell_dofs = acce_cells.at(cell_i).size();
                    const auto num_tot_dofs = lhs.rows();
                    const auto num_faces_dofs = num_tot_dofs - num_cell_dofs;

                    rhs.head(num_cell_dofs) = vector_type::Zero(num_cell_dofs);
                    lhs.topLeftCorner(num_cell_dofs, num_cell_dofs) =
                        matrix_type::Identity(num_cell_dofs, num_cell_dofs);
                    lhs.topRightCorner(num_faces_dofs, num_cell_dofs) *= 0.;
                    lhs.bottomLeftCorner(num_faces_dofs, num_cell_dofs) *= 0.;
                }
            }

            // std::cout << "R: " << rhs.norm() << std::endl;
            // std::cout << rhs.transpose() << std::endl;

            // Static Condensation
            // std::cout << "StatCond" << std::endl;
            tc.tic();
            const auto scnp = make_vector_static_condensation_withMatrix(msh, cl, degree_infos, lhs, rhs, true);

            m_AL[cell_i] = std::get<1>(scnp);
            m_bL[cell_i] = std::get<2>(scnp);

            tc.toc();
            ai.m_time_statcond += tc.elapsed();

            const auto &lc = std::get<0>(scnp);

            // std::cout << "rhs: " << lc.second.norm() << std::endl;
            // std::cout << lc.second.transpose() << std::endl;
            // std::cout << "lhs: " << lc.first.norm() << std::endl;

            // std::cout << "Assemb" << std::endl;
            m_assembler.assemble_nonlinear(msh, cl, bnd, lc.first, lc.second, depl_faces);
        }

        m_F_int = sqrt(m_F_int);
        // std::cout << "F_int: " << m_F_int << std::endl;

        m_assembler.impose_neumann_boundary_conditions(msh, bnd);
        m_assembler.finalize();

        ttot.toc();
        ai.m_time_assembly = ttot.elapsed();
        ai.m_linear_system_size = m_assembler.LHS.rows();
        return ai;
    }

    SolveInfo
    solve(void)
    {
        timecounter tc;

        tc.tic();
        m_system_displ = vector_type::Zero(m_assembler.LHS.rows());

#ifdef HAVE_INTEL_MKL
        solvers::pardiso_params<scalar_type> pparams;
        mkl_pardiso(pparams, m_assembler.LHS, m_assembler.RHS, m_system_displ);
#elif HAVE_MUMPS
        m_system_displ = mumps_lu(m_assembler.LHS, m_assembler.RHS);
#else
        throw std::runtime_error("No linear solver avalaible.");
#endif
        tc.toc();

        // std::cout << "LHS" << m_assembler.LHS << std::endl;
        // std::cout << "RHS" << m_assembler.RHS << std::endl;

        return SolveInfo(m_assembler.LHS.rows(), m_assembler.LHS.nonZeros(), tc.elapsed());
    }

    scalar_type
    postprocess(const mesh_type &msh,
                const bnd_type &bnd,
                const param_type &rp,
                const MeshDegreeInfo<mesh_type> &degree_infos,
                MultiTimeField<scalar_type> &fields)
    {
        timecounter tc;
        tc.tic();

        auto depl_faces = fields.getCurrentField(FieldName::DEPL_FACES);
        auto depl_cells = fields.getCurrentField(FieldName::DEPL_CELLS);
        auto depl = fields.getCurrentField(FieldName::DEPL);

        // Update cell
        for (auto &cl : msh)
        {
            const auto cell_i = msh.lookup(cl);

            const vector_type xdT =
                m_assembler.take_local_solution_nonlinear(msh, cl, bnd, m_system_displ, depl_faces);

            // static decondensation
            const vector_type xT = m_bL[cell_i] - m_AL[cell_i] * xdT;

            // Update element U^{i+1} = U^i + delta U^i
            depl.at(cell_i).head(xT.size()) += xT;
            depl.at(cell_i).tail(xdT.size()) += xdT;
            depl_cells.at(cell_i) += xT;

            // std::cout << "KT_F " << m_AL[cell_i].norm() << std::endl;
            // std::cout << "sol_F" << std::endl;
            // std::cout << xdT.transpose() << std::endl;
            // std::cout << "ft" << std::endl;
            // std::cout << m_bL[cell_i].transpose() << std::endl;
            // std::cout << "sol_T" << std::endl;
            // std::cout << xT.transpose() << std::endl;
            // std::cout << depl.at(cell_i).transpose() << std::endl;
        }
        fields.setCurrentField(FieldName::DEPL, depl);
        fields.setCurrentField(FieldName::DEPL_CELLS, depl_cells);

        // Update  unknowns
        // Update face Uf^{i+1} = Uf^i + delta Uf^i

        auto depl_faces_new = fields.getCurrentField(FieldName::DEPL_FACES);
        int face_i = 0;
        for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
        {
            const auto fc = *itor;
            depl_faces_new.at(face_i++) +=
                m_assembler.take_local_solution_nonlinear(msh, fc, bnd, m_system_displ, depl_faces);
        }
        fields.setCurrentField(FieldName::DEPL_FACES, depl_faces_new);

        m_dyna.postprocess(msh, m_time_step, fields);

        tc.toc();
        return tc.elapsed();
    }

    bool
    convergence(const param_type &rp, const size_t iter)
    {
        // norm of the solution
        auto error_un = std::sqrt(m_sol_norm);

        if (error_un <= scalar_type(10E-15))
        {
            error_un = scalar_type(10E16);
        }

        // norm of the rhs
        const scalar_type residual = m_assembler.RHS.norm();
        scalar_type max_error = 0.0;
        for (size_t i = 0; i < m_assembler.RHS.size(); i++)
            max_error = std::max(max_error, std::abs(m_assembler.RHS(i)));

        // norm of the increment
        const scalar_type error_incr = m_system_displ.norm();
        scalar_type relative_displ = error_incr / error_un;
        scalar_type relative_error = residual / m_F_int;

        if (iter == 0)
        {
            relative_displ = 1;
            relative_error = 1;
        }

        if (m_verbose)
        {
            std::string s_iter = "   " + std::to_string(iter) + "               ";
            s_iter.resize(9);

            if (iter == 0)
            {
                std::cout << "----------------------------------------------------------------------"
                             "------------------------"
                          << std::endl;
                std::cout << "| Iteration | Norme l2 incr | Relative incr |  Residual l2  | "
                             "Relative error | Maximum error |"
                          << std::endl;
                std::cout << "----------------------------------------------------------------------"
                             "------------------------"
                          << std::endl;
            }
            std::ios::fmtflags f(std::cout.flags());
            std::cout.precision(5);
            std::cout.setf(std::iostream::scientific, std::iostream::floatfield);
            std::cout << "| " << s_iter << " |   " << error_incr << " |   " << relative_displ << " |   " << residual
                      << " |   " << relative_error << "  |  " << max_error << "  |" << std::endl;
            std::cout << "-------------------------------------------------------------------------"
                         "---------------------"
                      << std::endl;
            std::cout.flags(f);
        }

        const scalar_type error = std::max(relative_displ, relative_error);

        if (!isfinite(error))
            throw std::runtime_error("Norm of residual is not finite");

        if (error <= rp.m_epsilon)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
};
}
} // end diskpp