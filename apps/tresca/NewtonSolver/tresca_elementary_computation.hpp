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

template<typename MeshType>
class tresca
{
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef disk::MaterialData<scalar_type>     material_type;
    typedef typename disk::hho_degree_info      hdi_type;

    const static int dimension = mesh_type::dimension;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    const mesh_type&     m_msh;
    const hdi_type&      m_hdi;
    const material_type& m_material_data;

  public:
    matrix_type K_int;
    vector_type RTF;
    vector_type F_int;

    tresca(const mesh_type& msh, const hdi_type& hdi, const material_type& material_data) :
      m_msh(msh), m_hdi(hdi), m_material_data(material_data)
    {
    }

    template<typename Function>
    void
    compute(const cell_type& cl, const Function& load, const matrix_type& ET, const vector_type& uTF)
    {
        const auto cell_degree = m_hdi.cell_degree();
        const auto grad_degree = m_hdi.grad_degree();
        const auto face_degree = m_hdi.face_degree();

        const auto cell_basis_size = disk::vector_basis_size(cell_degree, dimension, dimension);
        const auto grad_basis_size = disk::matrix_basis_size(grad_degree, dimension, dimension);
        const auto face_basis_size = disk::vector_basis_size(face_degree, dimension - 1, dimension);

        timecounter tc;

        const auto fcs            = faces(m_msh, cl);
        const auto num_faces      = fcs.size();
        const auto num_total_dofs = cell_basis_size + num_faces * face_basis_size;
        const auto dim_dofs       = dimension * dimension;

        const auto cb = make_vector_monomial_basis(m_msh, cl, cell_degree);

        // matrix_type AT = matrix_type::Zero(grad_basis_size, grad_basis_size);
        // vector_type aT = vector_type::Zero(grad_basis_size);

        RTF   = vector_type::Zero(num_total_dofs);
        F_int = vector_type::Zero(num_total_dofs);

        //   std::cout << "sol" << std::endl;
        //   std::cout << uTF.transpose() << std::endl;

        // std::cout << "AT: " << AT.norm() << std::endl;
        // std::cout << AT << std::endl;
        // std::cout << "at: " << aT.norm() << std::endl;

        // add volumic term
        RTF.head(cell_basis_size) = make_rhs(m_msh, cl, cb, load, 2);

        K_int = 2.0 * m_material_data.getMu() * ET;
        F_int = 2.0 * m_material_data.getMu() * ET * uTF;
        RTF -= F_int;

        //  std::cout << "K: " << K_int.norm() << std::endl;
        // // std::cout << K_int << std::endl;
        //  std::cout << "F: " << F_int.norm() << std::endl;

        assert(K_int.rows() == num_total_dofs);
        assert(K_int.cols() == num_total_dofs);
        assert(RTF.rows() == num_total_dofs);
    }
};

} // end namespace NLE