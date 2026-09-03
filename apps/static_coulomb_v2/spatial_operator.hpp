/*
 * spatial_operator.hpp : the u-independent, boundary-independent bulk.
 *
 * Owns the discretization (MeshDegreeInfo, LS: cell at k+1, faces at k) and
 * caches, per cell:
 *
 *     K[ci] = 2 mu gr.second + lam dr.second + 2 mu stab      (elasticity)
 *     G[ci] = gr.first                                        (gradrec)
 *     M[ci] = cell-block mass matrix                          (optional)
 *
 * Nitsche variants (hho_contact::trace_variant):
 *
 *   face : everything uniform, exactly the historical behaviour.
 *
 *   cell : Cascavita's cell version, in the realization of her reference
 *          Signorini Newton solver (solve_cells_full):
 *            - ALL face unknowns are KEPT (full uniform local layout);
 *            - on cells touching the contact boundary, the gradient and
 *              divergence reconstructions are built BLIND to the contact
 *              faces (their boundary integrals are skipped, their columns
 *              stay zero), so sigma_n does not depend on u_F there;
 *            - the HDG stabilization runs over ALL faces, which pins the
 *              contact-face unknown to the cell trace (no soft mode);
 *            - the Nitsche terms (contact_terms) read the CELL trace.
 *          Assembler, condensation and DOF layout are identical to the
 *          face version -- only K and G change, on the contact cells.
 */

#pragma once

#include <vector>
#include <functional>
#include <iostream>
#include <stdexcept>

#include "diskpp/methods/hho"
#include "diskpp/bases/bases.hpp"
#include "diskpp/adaptivity/adaptivity.hpp"

#include "hho_params.hpp"

namespace hho_contact {

template<typename Mesh>
class spatial_operator
{
public:
    using T           = typename Mesh::coordinate_type;
    using matrix_type = disk::dynamic_matrix<T>;
    using vector_type = disk::dynamic_vector<T>;
    using cell_type   = typename Mesh::cell_type;

    const Mesh&               msh;
    disk::MeshDegreeInfo<Mesh> di;

    size_t k;            // face degree
    size_t cell_deg;     // k+1 (Lehrenfeld-Schoberl)
    size_t cbs;          // cell basis size
    T      mu, lam;

    trace_variant variant;   // face (default) or cell Nitsche version

    std::vector<matrix_type> K;   // per-cell elasticity, full local dofs
    std::vector<matrix_type> G;   // per-cell gradient reconstruction
    std::vector<matrix_type> M;   // per-cell cell-block mass (empty if not requested)

    spatial_operator(const Mesh& m, size_t k_, T mu_, T lam_,
                     bool with_mass = false,
                     trace_variant variant_ = trace_variant::face,
                     std::function<bool(const typename Mesh::face_type&)>
                         is_contact_face = {})
        : msh(m), di(m, k_ + 1, k_), k(k_), cell_deg(k_ + 1),
          mu(mu_), lam(lam_), variant(variant_)
    {
        constexpr size_t DIM = Mesh::dimension;
        cbs = disk::vector_basis_size(cell_deg, DIM, DIM);
        const size_t fbs = disk::vector_basis_size(k, DIM - 1, DIM);

        if (variant == trace_variant::cell && !is_contact_face)
            throw std::invalid_argument(
                "spatial_operator: cell variant needs an is_contact_face "
                "predicate to know which faces the reconstruction must skip");

        const size_t nc = msh.cells_size();
        K.resize(nc);
        G.resize(nc);
        if (with_mass) M.resize(nc);

        size_t contact_cells = 0;
        size_t ci = 0;
        for (auto& cl : msh)
        {
            const auto fcs = faces(msh, cl);

            /* which of this cell's faces are on the contact boundary?      */
            std::vector<bool> skip(fcs.size(), false);
            bool has_contact = false;
            if (variant == trace_variant::cell)
                for (size_t i = 0; i < fcs.size(); i++)
                    if (is_contact_face(fcs[i]))
                    {
                        skip[i]     = true;
                        has_contact = true;
                    }

            if (!has_contact)
            {
                /* uniform path: byte-identical to the face version         */
                auto gr   = disk::make_matrix_hho_symmetric_gradrec(msh, cl, di);
                auto dr   = disk::make_hho_divergence_reconstruction(msh, cl, di);
                auto stab = disk::make_vector_hdg_stabilization(msh, cl, di);

                K[ci] = 2.0*mu*gr.second + lam*dr.second + 2.0*mu*stab;
                G[ci] = gr.first;

                {   // ¿la expansión conserva el operador? columnas de contacto han de ser 0
                    static bool once = false;
                    if (!once) {
                        once = true;
                        double ncontact = 0.0;
                        size_t full = cbs;
                        for (size_t i = 0; i < fcs.size(); i++) {
                            if (skip[i]) ncontact += G[ci].block(0, full, G[ci].rows(), fbs).norm();
                            full += fbs;
                        }
                        std::cout << "  [dbg-exp] G.rows=" << G[ci].rows()
                                << " G.cols=" << G[ci].cols()
                                << " nfull=" << (cbs + fcs.size()*fbs)
                                << " |G(:,contact)|=" << ncontact
                                << " |G|=" << G[ci].norm()
                                << " K.rows=" << K[ci].rows() << "\n";
                    }
                }
            }
            else
            {
                /* contact cell, cell variant: contact-blind reconstructions
                 * on the FULL layout.
                 *
                 * diskpp's builders skip a face when its DegreeInfo carries
                 * no unknowns, but then they also drop its columns.  So we
                 * call them with a per-cell CellDegreeInfo flagging the
                 * contact faces (compact, correct result) and pad the result
                 * back to the full layout with zero rows/cols at the contact
                 * face blocks -- identical to building them full-size with
                 * zero contact-face columns, i.e. to Karol's
                 * make_hho_contact_scalar_laplacian, in vectorial form.    */
                contact_cells++;

                std::vector<disk::DegreeInfo> finfos;
                finfos.reserve(fcs.size());
                for (size_t i = 0; i < fcs.size(); i++)
                    finfos.push_back(skip[i] ? disk::DegreeInfo(false, k)
                                             : disk::DegreeInfo(k));

                /* degrees mirror MeshDegreeInfo(msh, k+1, k):
                 * cell k+1, reconstruction k+1, gradient k                 */
                const disk::CellDegreeInfo<Mesh> cdi(disk::DegreeInfo(cell_deg),
                                                     disk::DegreeInfo(k + 1),
                                                     disk::DegreeInfo(k),
                                                     finfos);

                auto gr   = disk::make_matrix_hho_symmetric_gradrec(msh, cl, cdi);
                auto dr   = disk::make_hho_divergence_reconstruction(msh, cl, cdi);
                /* stabilization stays FULL: every face, contact included --
                 * this is what pins the contact-face unknown to the cell   */
                auto stab = disk::make_vector_hdg_stabilization(msh, cl, di);

                const size_t nfull = cbs + fcs.size()*fbs;
                const auto   map   = compact_to_full_map(skip, fbs);

                K[ci] = 2.0*mu*expand_square (gr.second, map, nfull)
                      + lam  *expand_square (dr.second, map, nfull)
                      + 2.0*mu*stab;
                G[ci] = expand_columns(gr.first, map, nfull);
            }

            if (with_mass)
            {
                auto cb = disk::make_vector_monomial_basis(msh, cl, cell_deg);
                M[ci]  = disk::make_mass_matrix(msh, cl, cb);
            }
            ci++;
        }

        if (variant == trace_variant::cell)
            std::cout << "  cell variant: contact-blind reconstruction on "
                      << contact_cells << " cells (face unknowns kept)\n";
    }

    size_t num_cells() const { return msh.cells_size(); }

    size_t faces_dofs(const cell_type& cl) const
    {
        return disk::vector_faces_dofs(msh, cl, di);
    }

    size_t total_dofs(const cell_type& cl) const
    {
        return cbs + faces_dofs(cl);
    }

private:
    /* compact index -> full index, cell block first, then faces in order   */
    std::vector<size_t>
    compact_to_full_map(const std::vector<bool>& skip, const size_t fbs) const
    {
        std::vector<size_t> map;
        map.reserve(cbs + skip.size()*fbs);

        for (size_t j = 0; j < cbs; j++)
            map.push_back(j);

        size_t full = cbs;
        for (size_t i = 0; i < skip.size(); i++)
        {
            if (!skip[i])
                for (size_t t = 0; t < fbs; t++)
                    map.push_back(full + t);
            full += fbs;
        }
        return map;
    }

    matrix_type
    expand_square(const matrix_type& Mc, const std::vector<size_t>& map,
                  const size_t nfull) const
    {
        matrix_type Mf = matrix_type::Zero(nfull, nfull);
        for (size_t i = 0; i < map.size(); i++)
            for (size_t j = 0; j < map.size(); j++)
                Mf(map[i], map[j]) = Mc(i, j);
        return Mf;
    }

    matrix_type
    expand_columns(const matrix_type& Gc, const std::vector<size_t>& map,
                   const size_t nfull) const
    {
        matrix_type Gf = matrix_type::Zero(Gc.rows(), nfull);
        for (size_t j = 0; j < map.size(); j++)
            Gf.col(map[j]) = Gc.col(j);
        return Gf;
    }
};

} // namespace hho_contact