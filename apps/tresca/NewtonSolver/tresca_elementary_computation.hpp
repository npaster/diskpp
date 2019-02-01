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
#include "mechanics/behaviors/laws/behaviorlaws.hpp"
#include "mechanics/deformation_tensors.hpp"
#include "methods/hho"
#include "quadratures/quadratures.hpp"

#include "timecounter.h"

namespace MK
{

template<typename T, int DIM, typename GradBasis>
static_matrix<T, DIM, DIM>
eval_stress(const dynamic_matrix<T>&     ET,
            const dynamic_vector<T>&     uTF,
            const GradBasis&             gb,
            const disk::MaterialData<T>& material_data,
            const point<T, DIM>          pt)
{
    const auto ET_coeff = ET * uTF;

    const auto gphi = gb.eval_functions(pt);

    const auto ET_u = disk::eval(ET_coeff, gphi);

    return 2.0 * material_data.getMu() * ET_u +
           material_data.getLambda() * ET_u.trace() * static_matrix<T, DIM, DIM>::Identity();
}

template<typename MeshType>
class tresca
{
  private:
    typedef MeshType                                   mesh_type;
    typedef typename mesh_type::coordinate_type        scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;
    typedef point<scalar_type, mesh_type::dimension>   point_type;
    typedef disk::MaterialData<scalar_type>            material_type;
    typedef typename disk::hho_degree_info             hdi_type;
    typedef ParamRun<scalar_type>                      param_type;
    typedef disk::BoundaryConditions<mesh_type, false> bnd_type;

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
    make_hho_aT(const cell_type& cl, const matrix_type& ET, const matrix_type& ST)
    {

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

            // be carefull (I:I)_F = d
            const auto qp_d            = qp.weight() * dimension;
            const auto qp_d_gphi_trace = disk::priv::inner_product(qp_d, gphi_trace);

            mass_div += disk::priv::outer_product(qp_d_gphi_trace, gphi_trace);
        }

        aT += m_material_data.getLambda() * ET.transpose() * mass_div * ET;

        // stabilization part
        if (m_rp.m_stab)
        {
            aT += m_rp.m_beta * ST;
        }

        return aT;
    }

    template<typename CellBasis>
    vector_type
    make_hho_uT_n(const vector_static& n, const CellBasis& cb, const point_type& pt)
    {
        const auto c_phi = cb.eval_functions(pt);

        //(phi_T . n)
        return disk::priv::inner_product(c_phi, n);
    }

    template<typename CellBasis>
    Matrix<scalar_type, Dynamic, dimension>
    make_hho_uT_t(const vector_static& n, const CellBasis& cb, const point_type& pt)
    {
        const auto c_phi = cb.eval_functions(pt);
        const auto uT_n  = make_hho_uT_n(n, cb, pt);

        // phi_T - (phi_T . n)n
        return c_phi - disk::priv::inner_product(disk::priv::inner_product(c_phi, n), n);
    }

    template<typename GradBasis>
    Matrix<scalar_type, Dynamic, dimension>
    make_hho_sigma_n(const matrix_type& ET, const vector_static& n, const GradBasis& gb, const point_type& pt)
    {
        Matrix<scalar_type, Dynamic, dimension> sigma_n =
          Matrix<scalar_type, Dynamic, dimension>::Zero(num_total_dofs, dimension);

        const auto gphi   = gb.eval_functions(pt);
        const auto gphi_n = disk::priv::inner_product(gphi, n);

        // std::cout << "g_phi_n: " << gphi_n.norm() << std::endl;
        // std::cout << gphi_n.transpose() << std::endl;

        sigma_n = 2.0 * m_material_data.getMu() * ET.transpose() * gphi_n;

        // std::cout << "sigma_n1: "  << std::endl;
        // std::cout << (2.0 * m_material_data.getMu() * ET.transpose() * gphi_n).transpose() << std::endl;

        const auto gphi_trace_n = disk::priv::inner_product(disk::trace(gphi), n);

        // std::cout << "gphi_trace_n: " << gphi_trace_n.norm() << std::endl;
        // std::cout << gphi_trace_n.transpose() << std::endl;

        sigma_n += m_material_data.getLambda() * ET.transpose() * gphi_trace_n;

        // std::cout << "sigma_n2: " << std::endl;
        // std::cout << (m_material_data.getLambda() * ET.transpose() * gphi_trace_n).transpose() << std::endl;

        // sigma. n
        return sigma_n;
    }

    template<typename GradBasis>
    vector_type
    make_hho_sigma_nn(const matrix_type& ET, const vector_static& n, const GradBasis& gb, const point_type& pt)
    {
        const auto sigma_n = make_hho_sigma_n(ET, n, gb, pt);

        // sigma_n . n
        return disk::priv::inner_product(sigma_n, n);
    }

    vector_type
    make_hho_sigma_nn(const Matrix<scalar_type, Dynamic, dimension>& sigma_n, const vector_static& n)
    {
        // sigma_n . n
        return disk::priv::inner_product(sigma_n, n);
    }

    template<typename GradBasis>
    Matrix<scalar_type, Dynamic, dimension>
    make_hho_sigma_nt(const matrix_type& ET, const vector_static& n, const GradBasis& gb, const point_type& pt)
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
                   const scalar_type& gamma_F)
    {
        vector_type phi_n = theta * sigma_nn;
        phi_n.head(cell_basis_size) -= gamma_F * uT_n;

        // theta * sigma_nn - gamma uT_n
        return phi_n;
    }

    Matrix<scalar_type, Dynamic, dimension>
    make_hho_phi_t(const matrix_type& sigma_nt,
                   const matrix_type& uT_t,
                   const scalar_type& theta,
                   const scalar_type& gamma_F)
    {
        auto phi_t = theta * sigma_nt;
        phi_t.block(0, 0, cell_basis_size, dimension) -= gamma_F * uT_t;

        // theta * sigma_nt - gamma uT_t
        return phi_t;
    }

    // compute theta/gamma *(sigma_n, sigma_n)_Fc
    matrix_type
    make_hho_nitsche(const cell_type& cl, const matrix_type& ET)
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

                    // std::cout << "sigma_n: " << sigma_n.norm() << std::endl;
                    // std::cout << sigma_n << std::endl;
                    // std::cout << "qp_sigma_n: " << qp_sigma_n.norm() << std::endl;
                    // std::cout << qp_sigma_n << std::endl;

                    nitsche += disk::priv::outer_product(qp_sigma_n, sigma_n);
                }
            }
        }

        return m_rp.m_theta * nitsche;
    }

    // compute (phi_n_theta, H(-phi_n_1(u))*phi_n_1)_FC / gamma
    matrix_type
    make_hho_heaviside_contact(const cell_type& cl, const matrix_type& ET, const vector_type& uTF)
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
                const auto qps     = integrate(m_msh, fc, 2 * qp_deg);
                const auto hF      = diameter(m_msh, fc);
                const auto gamma_F = m_rp.m_gamma_0 / hF;

                for (auto& qp : qps)
                {
                    const vector_type sigma_nn = make_hho_sigma_nn(ET, n, gb, qp.point());
                    const vector_type uT_n     = make_hho_uT_n(n, cb, qp.point());

                    const vector_type phi_n_1 = make_hho_phi_n(sigma_nn, uT_n, 1., gamma_F);

                    // Heaviside(-phi_n_1(u))
                    if (phi_n_1.dot(uTF) < scalar_type(0))
                    {
                        const vector_type phi_n_theta = make_hho_phi_n(sigma_nn, uT_n, m_rp.m_theta, gamma_F);
                        const auto qp_phi_n_theta     = disk::priv::inner_product(qp.weight() / gamma_F, phi_n_theta);

                        lhs += disk::priv::outer_product(qp_phi_n_theta, phi_n_1);
                    }
                }
            }
        }

        return lhs;
    }

    vector_type
    make_hho_negative_contact(const cell_type& cl, const matrix_type& ET, const vector_type& uTF)
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
                const auto qps     = integrate(m_msh, fc, 2 * qp_deg);
                const auto hF      = diameter(m_msh, fc);
                const auto gamma_F = m_rp.m_gamma_0 / hF;

                for (auto& qp : qps)
                {
                    const vector_type sigma_nn = make_hho_sigma_nn(ET, n, gb, qp.point());
                    const vector_type uT_n     = make_hho_uT_n(n, cb, qp.point());

                    const vector_type phi_n_1 = make_hho_phi_n(sigma_nn, uT_n, 1., gamma_F);

                    const scalar_type negative = phi_n_1.dot(uTF);
                    // [negative]_R-
                    if (negative < scalar_type(0))
                    {
                        const vector_type phi_n_theta = make_hho_phi_n(sigma_nn, uT_n, m_rp.m_theta, gamma_F);

                        rhs += (qp.weight() / gamma_F * negative) * phi_n_theta;
                    }
                }
            }
        }
        return rhs;
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

        const auto fcs       = faces(m_msh, cl);
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
};

} // end namespace NLE