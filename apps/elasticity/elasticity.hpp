/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#pragma once

#include "common/eigen.hpp"
#include "quadratures/quadratures.hpp"
#include "bases/bases.hpp"

#include "hho/hho.hpp"

namespace disk {


   template<typename Mesh, typename CellBasisType, typename CellQuadType,
    typename FaceBasisType, typename FaceQuadType, typename Solution>
class l2_error
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef CellBasisType                       cell_basis_type;
    typedef CellQuadType                        cell_quadrature_type;
    typedef FaceBasisType                       face_basis_type;
    typedef FaceQuadType                        face_quadrature_type;


    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;




    cell_basis_type                             cell_basis;
    cell_quadrature_type                        cell_quadrature;

    face_basis_type                             face_basis;
    face_quadrature_type                        face_quadrature;



    size_t                                      m_degree;

public:
    matrix_type     oper;
    matrix_type     data;

    l2_error()
        : m_degree(1)
    {
        cell_basis          = cell_basis_type(m_degree);
        cell_quadrature     = cell_quadrature_type(m_degree);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    l2_error(size_t degree)
        : m_degree(degree)
    {
        cell_basis          = cell_basis_type(m_degree);
        cell_quadrature     = cell_quadrature_type(m_degree);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }



    scalar_type compute(const mesh_type& msh, const cell_type& cl,
                 const vector_type& x, Solution sol)
    {

         disk::projector_elas<mesh_type,
                                                   cell_basis_type,
                                                   cell_quadrature_type,
                                                   face_basis_type,
                                                   face_quadrature_type> proj(m_degree);


            dynamic_vector<scalar_type> true_dof = proj.compute_cell(msh, cl, sol);
            dynamic_vector<scalar_type> comp_dof = x.block(0,0,true_dof.size(), 1);

            std::cout << " " << '\n';
            for (size_t i = 0; i < comp_dof.size(); i++) {
               std::cout << true_dof(i) <<" "<< comp_dof(i)  << '\n';
            }
            dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
            scalar_type err_dof = diff_dof.dot(proj.cell_mm * diff_dof);


            return sqrt(err_dof);

    }
};

template<typename MeshType>
void
test_cell_projector(const MeshType& msh)
{
   typedef MeshType                            mesh_type;

   typedef typename mesh_type::scalar_type     scalar_type;

   typedef typename mesh_type::cell            cell_type;
   typedef typename mesh_type::face            face_type;

   typedef static_vector<scalar_type, 2> result_type;

   typedef disk::quadrature<mesh_type, cell_type>   cell_quadrature_type;
   typedef disk::quadrature<mesh_type, face_type>   face_quadrature_type;
   typedef disk::scaled_monomial_vector_basis<mesh_type, cell_type>  cell_basis_type;
   typedef disk::scaled_monomial_vector_basis<mesh_type, face_type>  face_basis_type;
   typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>  divc_basis_type;
   typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>  divf_basis_type;

   const size_t degree = 3;

   disk::projector_elas<mesh_type,
                        cell_basis_type,
                        cell_quadrature_type,
                        face_basis_type,
                        face_quadrature_type> proj(degree);

   disk::projector<mesh_type,
                        divc_basis_type,
                        cell_quadrature_type,
                        divf_basis_type,
                        face_quadrature_type> proj2(degree);

   const size_t dim = msh.dimension;

   std::cout << "Test projector degree 3 dim 2" << '\n';

   auto sol = [](const point<scalar_type,2>& p) -> result_type {
      scalar_type fx = 0.0;
      scalar_type fy = 0.0;
      const size_t degree = 3;
      for (size_t i = 0; i <= degree; i++) {
         fx += (i+1) * iexp_pow(p.x(), i);
         fy += (i+2) * iexp_pow(p.y(), i);
      }
      return result_type{fx,fy};
   };


   auto sol2 = [](const point<scalar_type,2>& p) -> auto {
      scalar_type fx = 0.0;
      scalar_type fy = 0.0;
      const size_t degree = 3;
      for (size_t i = 0; i <= degree; i++) {
         fx += iexp_pow(p.x(), i);
      }
      return fx;
   };

   Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> true_dof2;
   true_dof2.resize( 10, 1 );

   for (size_t i = 0; i < 10; i++) {
      true_dof2(i) = 0.0;
   }

   true_dof2(0) = 1.0;
   true_dof2(1) = 1.0;
   true_dof2(3) = 1.0;
   true_dof2(6) = 1.0;

   for (auto& cl : msh)
   {
      dynamic_vector<scalar_type> comp_dof = proj2.compute_cell(msh, cl, sol2);
      auto diff_dof = comp_dof - true_dof2;
      std::cout << "comp_dof" << '\n';
      std::cout << comp_dof << '\n';
      std::cout << "ecart" << '\n';
      std::cout << diff_dof.dot(proj2.cell_mm*diff_dof) << '\n';
   }
}

} // namespace disk
