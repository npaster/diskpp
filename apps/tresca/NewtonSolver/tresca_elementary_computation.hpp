/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
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
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

#pragma once

#include <cassert>

#include "bases/bases.hpp"
#include "common/eigen.hpp"
#include "mechanics/behaviors/maths_tensor.hpp"
#include "methods/hho"
#include "quadratures/quadratures.hpp"

#include "timecounter.h"

namespace MK
{

template<typename Mesh>
std::pair<Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
          Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>>
make_matrix_symmetric_gradrec(const Mesh&                                   msh,
                              const typename Mesh::cell_type&               cl,
                              const disk::hho_degree_info&                  di,
                              const disk::vector_boundary_conditions<Mesh>& bnd)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const size_t N = Mesh::dimension;

    const auto graddeg = di.grad_degree();
    const auto celdeg  = di.cell_degree();
    const auto facdeg  = di.face_degree();

    const auto gb = make_sym_matrix_monomial_basis(msh, cl, graddeg);
    const auto cb = make_vector_monomial_basis(msh, cl, celdeg);

    const auto gbs = disk::sym_matrix_basis_size(graddeg, Mesh::dimension, Mesh::dimension);
    const auto cbs = disk::vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    const auto fbs = disk::vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension);

    const auto fcs       = bnd.faces_without_contact(cl);
    const auto num_faces = fcs.size();

    matrix_type gr_lhs = matrix_type::Zero(gbs, gbs);
    matrix_type gr_rhs = matrix_type::Zero(gbs, cbs + num_faces * fbs);

    // this is very costly to build it
    const auto qps = integrate(msh, cl, 2 * graddeg);

    size_t dec = 0;
    if (N == 3)
        dec = 6;
    else if (N == 2)
        dec = 3;
    else
        std::logic_error("Expected 3 >= dim > 1");

    for (auto& qp : qps)
    {
        const auto gphi = gb.eval_functions(qp.point());

        for (size_t j = 0; j < gbs; j++)
        {
            const auto qp_gphi_j = disk::priv::inner_product(qp.weight(), gphi[j]);
            for (size_t i = j; i < gbs; i += dec)
                gr_lhs(i, j) += disk::priv::inner_product(gphi[i], qp_gphi_j);
        }
    }

    // upper part
    for (size_t j = 0; j < gbs; j++)
        for (size_t i = 0; i < j; i++)
            gr_lhs(i, j) = gr_lhs(j, i);

    // compute rhs
    if (celdeg > 0)
    {
        const auto qpc = integrate(msh, cl, graddeg + celdeg - 1);
        for (auto& qp : qpc)
        {
            const auto gphi    = gb.eval_functions(qp.point());
            const auto dphi    = cb.eval_sgradients(qp.point());
            const auto qp_dphi = disk::priv::inner_product(qp.weight(), dphi);

            gr_rhs.block(0, 0, gbs, cbs) += disk::priv::outer_product(gphi, qp_dphi);

        } // end qp
    }

    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);
        const auto fb = make_vector_monomial_basis(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, graddeg + std::max(celdeg, facdeg));
        for (auto& qp : qps_f)
        {
            const auto cphi = cb.eval_functions(qp.point());
            const auto gphi = gb.eval_functions(qp.point());
            const auto fphi = fb.eval_functions(qp.point());

            const auto qp_gphi_n = disk::priv::inner_product(gphi, disk::priv::inner_product(qp.weight(), n));
            gr_rhs.block(0, cbs + i * fbs, gbs, fbs) += disk::priv::outer_product(qp_gphi_n, fphi);
            gr_rhs.block(0, 0, gbs, cbs) -= disk::priv::outer_product(qp_gphi_n, cphi);
        }
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_scalar_hdg_stabilization(const Mesh&                                   msh,
                              const typename Mesh::cell_type&               cl,
                              const disk::hho_degree_info&                  di,
                              const disk::vector_boundary_conditions<Mesh>& bnd)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();

    const auto cbs = disk::scalar_basis_size(celdeg, Mesh::dimension);
    const auto fbs = disk::scalar_basis_size(facdeg, Mesh::dimension - 1);

    const auto fcs        = bnd.faces_without_contact(cl);
    const auto num_faces  = fcs.size();
    const auto total_dofs = cbs + num_faces * fbs;

    matrix_type       data = matrix_type::Zero(total_dofs, total_dofs);
    const matrix_type If   = matrix_type::Identity(fbs, fbs);

    auto cb = make_scalar_monomial_basis(msh, cl, celdeg);

    for (size_t i = 0; i < num_faces; i++)
    {
        const auto fc = fcs[i];
        const auto h  = diameter(msh, fc);
        auto       fb = make_scalar_monomial_basis(msh, fc, facdeg);

        matrix_type oper  = matrix_type::Zero(fbs, total_dofs);
        matrix_type tr    = matrix_type::Zero(fbs, total_dofs);
        matrix_type mass  = make_mass_matrix(msh, fc, fb);
        matrix_type trace = matrix_type::Zero(fbs, cbs);

        oper.block(0, cbs + i * fbs, fbs, fbs) = -If;

        const auto qps = integrate(msh, fc, facdeg + celdeg);
        for (auto& qp : qps)
        {
            const auto c_phi = cb.eval_functions(qp.point());
            const auto f_phi = fb.eval_functions(qp.point());

            assert(c_phi.rows() == cbs);
            assert(f_phi.rows() == fbs);
            assert(c_phi.cols() == f_phi.cols());

            trace += disk::priv::outer_product(disk::priv::inner_product(qp.weight(), f_phi), c_phi);
        }

        tr.block(0, cbs + i * fbs, fbs, fbs) = -mass;
        tr.block(0, 0, fbs, cbs)             = trace;

        oper.block(0, 0, fbs, cbs) = mass.ldlt().solve(trace);
        data += oper.transpose() * tr * (1. / h);

        // std::cout << data << std::endl;
    }

    return data;
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_vector_hdg_stabilization(const Mesh&                                   msh,
                              const typename Mesh::cell_type&               cl,
                              const disk::hho_degree_info&                  di,
                              const disk::vector_boundary_conditions<Mesh>& bnd)
{
    const auto hdg_scalar_stab = make_scalar_hdg_stabilization(msh, cl, di, bnd);

    return disk::priv::compute_lhs_vector(msh, cl, di, hdg_scalar_stab, bnd);
}

template<typename Mesh, typename T>
auto
make_vector_static_condensation_withMatrix(const Mesh&                                                      msh,
                                           const typename Mesh::cell_type&                                  cl,
                                           const disk::hho_degree_info&                                     di,
                                           const disk::vector_boundary_conditions<Mesh>&                    bnd,
                                           const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                                           const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    const auto num_cell_dofs = disk::vector_basis_size(di.cell_degree(), Mesh::dimension, Mesh::dimension);
    const auto num_face_dofs = disk::vector_basis_size(di.face_degree(), Mesh::dimension - 1, Mesh::dimension);

    const auto fcs            = bnd.faces_without_contact(cl);
    const auto num_faces      = fcs.size();
    const auto num_faces_dofs = num_faces * num_face_dofs;
    const auto num_total_dofs = num_cell_dofs + num_faces_dofs;

    assert(lhs.rows() == lhs.cols());
    assert(lhs.cols() == num_total_dofs);

    if ((rhs.size() != num_cell_dofs) && (rhs.size() != num_total_dofs))
    {
        throw std::invalid_argument("static condensation: incorrect size of the rhs");
    }

    const matrix_type K_TT = lhs.topLeftCorner(num_cell_dofs, num_cell_dofs);
    const matrix_type K_TF = lhs.topRightCorner(num_cell_dofs, num_faces_dofs);
    const matrix_type K_FT = lhs.bottomLeftCorner(num_faces_dofs, num_cell_dofs);
    const matrix_type K_FF = lhs.bottomRightCorner(num_faces_dofs, num_faces_dofs);

    assert(K_TT.cols() == num_cell_dofs);
    assert(K_TT.cols() + K_TF.cols() == lhs.cols());
    assert(K_TT.rows() + K_FT.rows() == lhs.rows());
    assert(K_TF.rows() + K_FF.rows() == lhs.rows());
    assert(K_FT.cols() + K_FF.cols() == lhs.cols());

    const vector_type cell_rhs  = rhs.head(num_cell_dofs);
    vector_type       faces_rhs = vector_type::Zero(num_faces_dofs);

    if (rhs.size() == num_total_dofs)
    {
        faces_rhs = rhs.tail(num_faces_dofs);
    }

    const auto K_TT_lu = K_TT.lu();
    // if (K_TT_lu.info() != Eigen::Success)
    // {
    //     throw std::invalid_argument("static condensation: K_TT is not positive definite");
    // }

    const matrix_type AL = K_TT_lu.solve(K_TF);
    const vector_type bL = K_TT_lu.solve(cell_rhs);

    const matrix_type AC = K_FF - K_FT * AL;
    const vector_type bC = faces_rhs - K_FT * bL;

    return std::make_tuple(std::make_pair(AC, bC), AL, bL);
}

template<typename MeshType>
class tresca
{
  private:
    typedef MeshType                                    mesh_type;
    typedef typename mesh_type::coordinate_type         scalar_type;
    typedef typename mesh_type::cell                    cell_type;
    typedef typename mesh_type::face                    face_type;
    typedef point<scalar_type, mesh_type::dimension>    point_type;
    typedef disk::MaterialData<scalar_type>             material_type;
    typedef typename disk::hho_degree_info              hdi_type;
    typedef ParamRun<scalar_type>                       param_type;
    typedef disk::vector_boundary_conditions<mesh_type> bnd_type;

    const static int dimension = mesh_type::dimension;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    typedef static_matrix<scalar_type, dimension, dimension> matrix_static;
    typedef static_vector<scalar_type, dimension>            vector_static;

    const mesh_type&     m_msh;
    const hdi_type&      m_hdi;
    const material_type& m_material_data;
    const param_type&    m_rp;
    const bnd_type&      m_bnd;

    size_t cell_basis_size, grad_basis_size, face_basis_size, num_total_dofs;

    // contact contrib;

    matrix_type
    make_hho_aT(const cell_type& cl, const matrix_type& ET, const matrix_type& ST) const
    {
        // grad_sym part
        matrix_type mass_grad = matrix_type::Zero(grad_basis_size, grad_basis_size);
        const auto  gb        = make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());

        // this is very costly to build it
        const auto qps = integrate(m_msh, cl, 2 * m_hdi.grad_degree());

        size_t dec = 0;
        if (dimension == 3)
            dec = 6;
        else if (dimension == 2)
            dec = 3;
        else
            std::logic_error("Expected 3 >= dim > 1");

        for (auto& qp : qps)
        {
            const auto gphi = gb.eval_functions(qp.point());

            for (size_t j = 0; j < grad_basis_size; j++)
            {
                const auto qp_gphi_j = disk::priv::inner_product(qp.weight(), gphi[j]);
                for (size_t i = j; i < grad_basis_size; i += dec)
                    mass_grad(i, j) += disk::priv::inner_product(gphi[i], qp_gphi_j);
            }
        }

        // upper part
        for (size_t j = 0; j < grad_basis_size; j++)
            for (size_t i = 0; i < j; i++)
                mass_grad(i, j) = mass_grad(j, i);

        matrix_type aT = 2.0 * m_material_data.getMu() * ET.transpose() * mass_grad * ET;

        // divergence part
        matrix_type mass_div = matrix_type::Zero(grad_basis_size, grad_basis_size);

        for (auto& qp : qps)
        {
            const auto gphi       = gb.eval_functions(qp.point());
            const auto gphi_trace = disk::trace(gphi);

            const auto qp_gphi_trace = disk::priv::inner_product(qp.weight(), gphi_trace);

            mass_div += disk::priv::outer_product(qp_gphi_trace, gphi_trace);
        }

        aT += m_material_data.getLambda() * ET.transpose() * mass_div * ET;

        // stabilization part
        if (m_rp.m_stab)
        {
            aT += m_rp.m_beta * ST;
        }

        return aT;
    }

    scalar_type
    make_hho_distance(const point_type& pt, const vector_static& n) const
    {
        if (std::abs(n(0)) < scalar_type(1E-12))
            return pt.y();

        scalar_type delta = n(1) / n(0);
        if (std::abs(delta) < scalar_type(1E-12))
            return 10E12;

        return pt.y() * (delta * n(0) + scalar_type(1) * n(1));
    }

    template<typename CellBasis>
    vector_type
    make_hho_uT_n(const vector_static& n, const CellBasis& cb, const point_type& pt) const
    {
        const auto c_phi = cb.eval_functions(pt);

        //(phi_T . n)
        return disk::priv::inner_product(c_phi, n);
    }

    template<typename CellBasis>
    Matrix<scalar_type, Dynamic, dimension>
    make_hho_uT_t(const vector_static& n, const CellBasis& cb, const point_type& pt) const
    {
        const auto c_phi = cb.eval_functions(pt);
        const auto uT_n  = make_hho_uT_n(n, cb, pt);

        // phi_T - (phi_T . n)n
        return c_phi - disk::priv::inner_product(disk::priv::inner_product(c_phi, n), n);
    }

    template<typename GradBasis>
    Matrix<scalar_type, Dynamic, dimension>
    make_hho_sigma_n(const matrix_type& ET, const vector_static& n, const GradBasis& gb, const point_type& pt) const
    {
        Matrix<scalar_type, Dynamic, dimension> sigma_n =
          Matrix<scalar_type, Dynamic, dimension>::Zero(num_total_dofs, dimension);

        const auto gphi   = gb.eval_functions(pt);
        const auto gphi_n = disk::priv::inner_product(gphi, n);

        sigma_n = 2.0 * m_material_data.getMu() * ET.transpose() * gphi_n;

        const auto gphi_trace_n = disk::priv::inner_product(disk::trace(gphi), n);

        sigma_n += m_material_data.getLambda() * ET.transpose() * gphi_trace_n;

        // sigma. n
        return sigma_n;
    }

    template<typename GradBasis>
    vector_type
    make_hho_sigma_nn(const matrix_type& ET, const vector_static& n, const GradBasis& gb, const point_type& pt) const
    {
        const auto sigma_n = make_hho_sigma_n(ET, n, gb, pt);

        // sigma_n . n
        return disk::priv::inner_product(sigma_n, n);
    }

    vector_type
    make_hho_sigma_nn(const Matrix<scalar_type, Dynamic, dimension>& sigma_n, const vector_static& n) const
    {
        // sigma_n . n
        return disk::priv::inner_product(sigma_n, n);
    }

    template<typename GradBasis>
    Matrix<scalar_type, Dynamic, dimension>
    make_hho_sigma_nt(const matrix_type& ET, const vector_static& n, const GradBasis& gb, const point_type& pt) const
    {
        const auto sigma_n  = make_hho_sigma_n(ET, n, gb, pt);
        const auto sigma_nn = make_hho_sigma_nn(sigma_n, n);

        // sigma_n - sigma_nn * n
        return sigma_n - disk::priv::inner_product(sigma_nn, n);
    }

    vector_type
    make_hho_phi_n(const vector_type& sigma_nn,
                   const vector_type& uT_n,
                   const scalar_type& theta,
                   const scalar_type& gamma_F) const
    {
        vector_type phi_n = theta * sigma_nn;
        phi_n.head(cell_basis_size) -= gamma_F * uT_n;

        // theta * sigma_nn - gamma uT_n
        return phi_n;
    }

    Matrix<scalar_type, Dynamic, dimension>
    make_hho_phi_t(const Matrix<scalar_type, Dynamic, dimension>& sigma_nt,
                   const Matrix<scalar_type, Dynamic, dimension>& uT_t,
                   const scalar_type&                             theta,
                   const scalar_type&                             gamma_F) const
    {
        Matrix<scalar_type, Dynamic, dimension> phi_t = theta * sigma_nt;
        phi_t.block(0, 0, cell_basis_size, dimension) -= gamma_F * uT_t;

        // theta * sigma_nt - gamma uT_t
        return phi_t;
    }

    // projection on the ball of radius alpha centered on 0
    vector_static
    make_proj_alpha(const vector_static& x, const scalar_type& alpha) const
    {
        const scalar_type x_norm = x.norm();

        if (x_norm <= alpha)
        {
            return x;
        }

        return alpha * x / x_norm;
    }

    // derivative of the projection on the ball of radius alpha centered on 0
    matrix_static
    make_d_proj_alpha(const vector_static& x, const scalar_type& alpha) const
    {
        const scalar_type x_norm = x.norm();

        if (x_norm <= alpha)
        {
            return matrix_static::Identity();
        }

        return alpha / x_norm * (matrix_static::Identity() - disk::Kronecker(x, x) / (x_norm * x_norm));
    }

    // compute theta/gamma *(sigma_n, sigma_n)_Fc
    matrix_type
    make_hho_nitsche(const cell_type& cl, const matrix_type& ET) const
    {
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());

        matrix_type nitsche = matrix_type::Zero(num_total_dofs, num_total_dofs);

        const auto fcs = faces(m_msh, cl);
        for (auto& fc : fcs)
        {
            if (m_bnd.is_contact_face(fc))
            {
                const auto n       = normal(m_msh, cl, fc);
                const auto qps     = integrate(m_msh, fc, 2 * m_hdi.grad_degree());
                const auto hF      = diameter(m_msh, fc);
                const auto gamma_F = m_rp.m_gamma_0 / hF;

                for (auto& qp : qps)
                {
                    const auto sigma_n    = make_hho_sigma_n(ET, n, gb, qp.point());
                    const auto qp_sigma_n = disk::priv::inner_product(qp.weight() / gamma_F, sigma_n);

                    nitsche += disk::priv::outer_product(qp_sigma_n, sigma_n);
                }
            }
        }

        return m_rp.m_theta * nitsche;
    }

    // compute (phi_n_theta, H(-phi_n_1(u))*phi_n_1)_FC / gamma
    matrix_type
    make_hho_heaviside_contact(const cell_type& cl, const matrix_type& ET, const vector_type& uTF) const
    {
        const auto cb = make_vector_monomial_basis(m_msh, cl, m_hdi.cell_degree());
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());

        matrix_type lhs = matrix_type::Zero(num_total_dofs, num_total_dofs);

        const auto fcs = faces(m_msh, cl);
        for (auto& fc : fcs)
        {
            if (m_bnd.is_contact_face(fc))
            {
                const auto n       = normal(m_msh, cl, fc);
                const auto qp_deg  = std::max(m_hdi.cell_degree(), m_hdi.grad_degree());
                const auto qps     = integrate(m_msh, fc, 2 * qp_deg + 2);
                const auto hF      = diameter(m_msh, fc);
                const auto gamma_F = m_rp.m_gamma_0 / hF;

                for (auto& qp : qps)
                {
                    const vector_type sigma_nn = make_hho_sigma_nn(ET, n, gb, qp.point());
                    const vector_type uT_n     = make_hho_uT_n(n, cb, qp.point());

                    const scalar_type phi_n_1_u = eval_phi_n(ET, gb, cb, uTF, n, gamma_F, qp.point());

                    // Heaviside(-phi_n_1(u))
                    if (phi_n_1_u <= scalar_type(0))
                    {
                        const vector_type phi_n_theta = make_hho_phi_n(sigma_nn, uT_n, m_rp.m_theta, gamma_F);
                        const vector_type phi_n_1     = make_hho_phi_n(sigma_nn, uT_n, scalar_type(1), gamma_F);
                        const auto qp_phi_n_theta     = disk::priv::inner_product(qp.weight() / gamma_F, phi_n_theta);

                        lhs += disk::priv::outer_product(qp_phi_n_theta, phi_n_1);
                    }
                }
            }
        }

        return lhs;
    }

    // compute (phi_n_theta, [phi_n_1(u)]R-)_FC / gamma
    vector_type
    make_hho_negative_contact(const cell_type& cl, const matrix_type& ET, const vector_type& uTF) const
    {
        const auto cb = make_vector_monomial_basis(m_msh, cl, m_hdi.cell_degree());
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());

        vector_type rhs = vector_type::Zero(num_total_dofs);

        const auto fcs = faces(m_msh, cl);
        for (auto& fc : fcs)
        {
            if (m_bnd.is_contact_face(fc))
            {
                const auto n       = normal(m_msh, cl, fc);
                const auto qp_deg  = std::max(m_hdi.cell_degree(), m_hdi.grad_degree());
                const auto qps     = integrate(m_msh, fc, 2 * qp_deg + 2);
                const auto hF      = diameter(m_msh, fc);
                const auto gamma_F = m_rp.m_gamma_0 / hF;

                for (auto& qp : qps)
                {
                    const vector_type sigma_nn = make_hho_sigma_nn(ET, n, gb, qp.point());
                    const vector_type uT_n     = make_hho_uT_n(n, cb, qp.point());

                    const scalar_type phi_n_1_u = eval_phi_n(ET, gb, cb, uTF, n, gamma_F, qp.point());

                    // [phi_n_1_u]_R-
                    if (phi_n_1_u <= scalar_type(0))
                    {
                        const vector_type phi_n_theta = make_hho_phi_n(sigma_nn, uT_n, m_rp.m_theta, gamma_F);

                        rhs += (qp.weight() / gamma_F * phi_n_1_u) * phi_n_theta;
                    }
                }
            }
        }
        return rhs;
    }

    // compute (phi_t_theta, [phi_t_1(u)]_(s))_FC / gamma
    vector_type
    make_hho_threshold_tresca(const cell_type& cl, const matrix_type& ET, const vector_type& uTF) const
    {
        const auto cb = make_vector_monomial_basis(m_msh, cl, m_hdi.cell_degree());
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());

        vector_type rhs = vector_type::Zero(num_total_dofs);

        const auto fcs = faces(m_msh, cl);
        for (auto& fc : fcs)
        {
            if (m_bnd.is_contact_face(fc))
            {
                const auto n       = normal(m_msh, cl, fc);
                const auto qp_deg  = std::max(m_hdi.cell_degree(), m_hdi.grad_degree());
                const auto qps     = integrate(m_msh, fc, 2 * qp_deg + 2);
                const auto hF      = diameter(m_msh, fc);
                const auto gamma_F = m_rp.m_gamma_0 / hF;

                for (auto& qp : qps)
                {
                    const auto sigma_nt = make_hho_sigma_nt(ET, n, gb, qp.point());
                    const auto uT_t     = make_hho_uT_t(n, cb, qp.point());

                    const auto phi_t_theta = make_hho_phi_t(sigma_nt, uT_t, m_rp.m_theta, gamma_F);

                    const vector_static phi_t_1_u_proj = eval_proj_phi_t(ET, gb, cb, uTF, n, gamma_F, qp.point());

                    const vector_static qp_phi_t_1_u_pro = qp.weight() * phi_t_1_u_proj / gamma_F;

                    rhs += disk::priv::inner_product(phi_t_theta, qp_phi_t_1_u_pro);
                }
            }
        }
        return rhs;
    }

    // compute (phi_t_theta, (d_proj_alpha(u)) phi_t_1)_FC / gamma
    matrix_type
    make_hho_matrix_tresca(const cell_type& cl, const matrix_type& ET, const vector_type& uTF) const
    {
        const auto cb = make_vector_monomial_basis(m_msh, cl, m_hdi.cell_degree());
        const auto gb = make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());

        matrix_type lhs = matrix_type::Zero(num_total_dofs, num_total_dofs);

        const auto fcs = faces(m_msh, cl);
        for (auto& fc : fcs)
        {
            if (m_bnd.is_contact_face(fc))
            {
                const auto n       = normal(m_msh, cl, fc);
                const auto qp_deg  = std::max(m_hdi.cell_degree(), m_hdi.grad_degree());
                const auto qps     = integrate(m_msh, fc, 2 * qp_deg + 2);
                const auto hF      = diameter(m_msh, fc);
                const auto gamma_F = m_rp.m_gamma_0 / hF;

                for (auto& qp : qps)
                {
                    const auto sigma_nt = make_hho_sigma_nt(ET, n, gb, qp.point());
                    const auto uT_t     = make_hho_uT_t(n, cb, qp.point());

                    const auto phi_t_1     = make_hho_phi_t(sigma_nt, uT_t, scalar_type(1), gamma_F);
                    const auto phi_t_theta = make_hho_phi_t(sigma_nt, uT_t, m_rp.m_theta, gamma_F);

                    const auto phi_t_1_u      = eval_phi_t(ET, gb, cb, uTF, n, gamma_F, qp.point());
                    const auto d_proj_phi_t_u = make_d_proj_alpha(phi_t_1_u, m_rp.m_threshold);

                    const auto d_proj_u_phi_t_1 = disk::priv::inner_product(d_proj_phi_t_u, phi_t_1);

                    const auto qp_phi_t_theta = disk::priv::inner_product(qp.weight() / gamma_F, phi_t_theta);

                    lhs += disk::priv::outer_product(qp_phi_t_theta, d_proj_u_phi_t_1);
                }
            }
        }

        return lhs;
    }

  public:
    matrix_type K_int;
    vector_type RTF;
    vector_type F_int;
    double      time_contact;

    tresca(const mesh_type&     msh,
           const hdi_type&      hdi,
           const material_type& material_data,
           const param_type&    rp,
           const bnd_type&      bnd) :
      m_msh(msh),
      m_hdi(hdi), m_material_data(material_data), m_rp(rp), m_bnd(bnd)
    {
        cell_basis_size = disk::vector_basis_size(m_hdi.cell_degree(), dimension, dimension);
        grad_basis_size = disk::sym_matrix_basis_size(m_hdi.grad_degree(), dimension, dimension);
        face_basis_size = disk::vector_basis_size(m_hdi.face_degree(), dimension - 1, dimension);
        num_total_dofs  = 0;
    }

    template<typename Function>
    void
    compute(const cell_type&   cl,
            const Function&    load,
            const matrix_type& ET,
            const matrix_type& ST,
            const vector_type& uTF,
            const bool&        has_vector_face)
    {
        timecounter tc;

        const auto fcs       = m_bnd.faces_without_contact(cl);
        const auto num_faces = fcs.size();
        num_total_dofs       = cell_basis_size + num_faces * face_basis_size;

        const auto cb = make_vector_monomial_basis(m_msh, cl, m_hdi.cell_degree());

        RTF   = vector_type::Zero(num_total_dofs);
        F_int = vector_type::Zero(num_total_dofs);
        K_int = matrix_type::Zero(num_total_dofs, num_total_dofs);

        // volumic contribution
        matrix_type AT = make_hho_aT(cl, ET, ST);

        // add volumic term
        RTF.head(cell_basis_size) += make_rhs(m_msh, cl, cb, load, 2);

        // contact contribution
        time_contact = 0.0;
        if (has_vector_face)
        {
            tc.tic();
            // compute theta/gamma *(sigma_n, sigma_n)_Fc
            AT -= make_hho_nitsche(cl, ET);

            // std::cout << "Nitche: " << std::endl;
            // std::cout << make_hho_nitsche(cl, ET) << std::endl;

            // compute (phi_n_theta, H(-phi_n_1(u))*phi_n_1)_FC / gamma
            K_int += make_hho_heaviside_contact(cl, ET, uTF);

            // std::cout << "Heaviside: " << std::endl;
            // std::cout << make_hho_heaviside_contact(cl, ET, uTF) << std::endl;

            // compute (phi_n_theta, [phi_n_1(u)]R-)_FC / gamma
            F_int += make_hho_negative_contact(cl, ET, uTF);

            // std::cout << "Negative: " << std::endl;
            // std::cout << make_hho_negative_contact(cl, ET, uTF).transpose() << std::endl;

            // friction contribution
            if (m_rp.m_frot)
            {
                // compute (phi_t_theta, [phi_t_1(u)]_s)_FC / gamma
                F_int += make_hho_threshold_tresca(cl, ET, uTF);

                // std::cout << "Threshold: " << std::endl;
                // std::cout << make_hho_threshold_tresca(cl, ET, uTF).transpose() << std::endl;

                // compute (phi_t_theta, (d_proj_alpha(u)) phi_t_1)_FC / gamma
                K_int += make_hho_matrix_tresca(cl, ET, uTF);
            }

            tc.toc();
            time_contact = tc.to_double();
        }

        K_int += AT;
        F_int += AT * uTF;
        RTF -= F_int;

        //  std::cout << "K: " << K_int.norm() << std::endl;
        //  std::cout << K_int << std::endl;
        //  std::cout << "R: " << RTF.norm() << std::endl;
        //  std::cout << RTF.transpose() << std::endl;

        assert(K_int.rows() == num_total_dofs);
        assert(K_int.cols() == num_total_dofs);
        assert(RTF.rows() == num_total_dofs);
    }

    template<typename CellBasis>
    scalar_type
    eval_uT_n(const CellBasis& cb, const vector_type& uTF, const vector_static& n, const point_type& pt) const
    {
        const vector_type uT_n = make_hho_uT_n(n, cb, pt);

        return uT_n.dot(uTF.head(cb.size()));
    }

    template<typename CellBasis>
    vector_static
    eval_uT_t(const CellBasis& cb, const vector_type& uTF, const vector_static& n, const point_type& pt) const
    {
        const auto uT_t = make_hho_uT_t(n, cb, pt);

        return uT_t.transpose() * (uTF.head(cb.size()));
    }

    template<typename GradBasis>
    scalar_type
    eval_sigma_nn(const matrix_type&   ET,
                  const GradBasis&     gb,
                  const vector_type&   uTF,
                  const vector_static& n,
                  const point_type&    pt) const
    {
        const vector_type sigma_nn = make_hho_sigma_nn(ET, n, gb, pt);

        return sigma_nn.dot(uTF);
    }

    template<typename GradBasis>
    vector_static
    eval_sigma_nt(const matrix_type&   ET,
                  const GradBasis&     gb,
                  const vector_type&   uTF,
                  const vector_static& n,
                  const point_type&    pt) const
    {
        const auto sigma_nt = make_hho_sigma_nt(ET, n, gb, pt);

        return sigma_nt.transpose() * uTF;
    }

    template<typename GradBasis, typename CellBasis>
    scalar_type
    eval_phi_n(const matrix_type&   ET,
               const GradBasis&     gb,
               const CellBasis&     cb,
               const vector_type&   uTF,
               const vector_static& n,
               const scalar_type&   gamma_F,
               const point_type&    pt) const
    {
        const vector_type sigma_nn = make_hho_sigma_nn(ET, n, gb, pt);
        const vector_type uT_n     = make_hho_uT_n(n, cb, pt);

        const vector_type phi_n_1   = make_hho_phi_n(sigma_nn, uT_n, scalar_type(1), gamma_F);
        const scalar_type phi_n_1_u = phi_n_1.dot(uTF);

        return phi_n_1_u;
    }

    template<typename GradBasis, typename CellBasis>
    scalar_type
    eval_proj_phi_n(const matrix_type&   ET,
                    const GradBasis&     gb,
                    const CellBasis&     cb,
                    const vector_type&   uTF,
                    const vector_static& n,
                    const scalar_type&   gamma_F,
                    const point_type&    pt) const
    {
        const scalar_type phi_n_1_u = eval_phi_n(ET, gb, cb, uTF, n, gamma_F, pt);

        if (phi_n_1_u <= scalar_type(0))
            return phi_n_1_u;

        return scalar_type(0);
    }

    template<typename GradBasis, typename CellBasis>
    vector_static
    eval_phi_t(const matrix_type&   ET,
               const GradBasis&     gb,
               const CellBasis&     cb,
               const vector_type&   uTF,
               const vector_static& n,
               const scalar_type&   gamma_F,
               const point_type&    pt) const
    {
        const auto sigma_nt = make_hho_sigma_nt(ET, n, gb, pt);
        const auto uT_t     = make_hho_uT_t(n, cb, pt);

        const auto          phi_t_1   = make_hho_phi_t(sigma_nt, uT_t, scalar_type(1), gamma_F);
        const vector_static phi_t_1_u = phi_t_1.transpose() * uTF;

        return phi_t_1_u;
    }

    template<typename GradBasis, typename CellBasis>
    vector_static
    eval_proj_phi_t(const matrix_type&   ET,
                    const GradBasis&     gb,
                    const CellBasis&     cb,
                    const vector_type&   uTF,
                    const vector_static& n,
                    const scalar_type&   gamma_F,
                    const point_type&    pt) const
    {
        const vector_static phi_t_1_u = eval_phi_t(ET, gb, cb, uTF, n, gamma_F, pt);

        return make_proj_alpha(phi_t_1_u, m_rp.m_threshold);
    }

    template<typename GradBasis>
    matrix_static
    eval_stress(const vector_type& ET_uTF, const GradBasis& gb, const point_type pt) const
    {
        const auto gphi = gb.eval_functions(pt);

        const auto ET_u = disk::eval(ET_uTF, gphi);

        return 2.0 * m_material_data.getMu() * ET_u +
               m_material_data.getLambda() * ET_u.trace() * matrix_static::Identity();
    }

    template<typename GradBasis>
    matrix_static
    eval_stress(const matrix_type& ET, const vector_type& uTF, const GradBasis& gb, const point_type pt) const
    {
        const vector_type ET_uTF = ET * uTF;

        return eval_stress(ET_uTF, gb, pt);
    }
};

} // end namespace MK