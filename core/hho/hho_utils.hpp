/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
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

#include "common/eigen.hpp"

namespace disk {

   namespace hho{

      template<typename T, int DIM>
      static_matrix<T,DIM, DIM>
      eval_gradient(const dynamic_vector<T>& gradrec_coeff,
                    const std::vector<static_matrix<T,DIM, DIM> >& base_grad)
      {
         static_matrix<T,DIM, DIM> ret = static_matrix<T,DIM, DIM>::Zero();
         assert(gradrec_coeff.size() == base_grad.size());

         for(size_t i = 0; i < base_grad.size(); i++){
            ret += gradrec_coeff(i) * base_grad[i];
         }

         return ret;
      }

      template<typename T, int DIM>
      static_vector<T,DIM>
      eval_gradient(const dynamic_vector<T>& gradrec_coeff,
                    const std::vector<static_vector<T,DIM> >& base_grad)
      {
         static_vector<T,DIM> ret = static_vector<T,DIM>::Zero();
         assert(gradrec_coeff.size() == base_grad.size());

         for(size_t i = 0; i < base_grad.size(); i++){
            ret += gradrec_coeff(i) * base_grad[i];
         }

         return ret;
      }



      template<typename T, int DIM>
      static_matrix<T,DIM, DIM>
      eval(const dynamic_vector<T>& tab_coeff,
            const std::vector<static_matrix<T,DIM, DIM> >& base)
      {
         static_matrix<T,DIM, DIM> ret = static_matrix<T,DIM, DIM>::Zero();
         assert(tab_coeff.size() == base.size());

         for(size_t i = 0; i < base.size(); i++){
            ret += tab_coeff(i) * base[i];
         }

         return ret;
      }

      template<typename T, int DIM>
      static_vector<T,DIM>
      eval(const dynamic_vector<T>& tab_coeff,
                    const std::vector<static_vector<T,DIM> >& base)
      {
         static_vector<T,DIM> ret = static_vector<T,DIM>::Zero();
         assert(tab_coeff.size() == base.size());

         for(size_t i = 0; i < base.size(); i++){
            ret += tab_coeff(i) * base[i];
         }

         return ret;
      }


/*
      //compute mass matrix
      namespace priv {
         template< typename BQData, typename BasisType>
         struct mass_matrix_F {
            typedef typename BQData::mesh_type                    mesh_type;
            typedef typename mesh_type::scalar_type               scalar_type;
            typedef typename mesh_type::cell                      cell_type;
            typedef typename mesh_type::face                      face_type;

            static void impl(const mesh_type& msh, const cell_type& cl, const BQData& bqd)
            { static_assert(sizeof(BQData) = -1  && "BQData not known in mass_matrix_F"); }

            static void impl(const mesh_type& msh, const face_type& fc, const BQData& bqd)
            { static_assert(sizeof(BQData) = -1  && "BQData not known in mass_matrix_F"); }
         };

         //scaled_monomial_scalar_basis
         template<typename BQData>
         struct mass_matrix_F<BQData, scaled_monomial_scalar_basis> {

            typedef typename BQData::mesh_type                    mesh_type;
            typedef typename mesh_type::scalar_type               scalar_type;
            typedef typename mesh_type::cell                      cell_type;
            typedef typename mesh_type::cell                      face_type;

            typedef dynamic_matrix<scalar_type>                   matrix_type;

            static matrix_type impl(const mesh_type& msh, const cell_type& cl, const BQData& bqd)
            {
               const auto cell_degree = bqd.cell_degree();
               const auto num_cell_dofs = howmany_dofs(bqd.cell_basis);

               matrix_type mass = matrix_type::Zero(num_cell_dofs, num_cell_dofs);

               const auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
               for (auto& qp : cell_quadpoints)
               {
               const matrix_type cphi = bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
                  assert(cphi.rows() == num_cell_dofs);

                  mass += qp.weight() * cphi * cphi.transpose();
               }

               return mass;
            }

            static matrix_type impl(const mesh_type& msh, const face_type& fc, const BQData& bqd)
            {
               const auto face_degree = bqd.face_degree();
               const auto num_face_dofs = howmany_dofs(bqd.face_basis);

               matrix_type mass = matrix_type::Zero(num_face_dofs, num_face_dofs);

               const auto face_quadpoints = bqd.face_quadrature.integrate(msh, fc);
               for (auto& qp : face_quadpoints)
               {
                  const matrix_type fphi = bqd.face_basis.eval_functions(msh, fc, qp.point(), 0, face_degree);
                  assert(fphi.rows() == num_face_dofs);

                  mass += qp.weight() * fphi * fphi.transpose();
               }

               return mass;
            }
         };

         //scaled_monomial_vector_basis
         template<typename BQData>
         struct mass_matrix_F<BQData, scaled_monomial_vector_basis> {

            typedef typename BQData::mesh_type                    mesh_type;
            typedef typename mesh_type::scalar_type               scalar_type;
            typedef typename mesh_type::cell                      cell_type;
            typedef typename mesh_type::face                      face_type

            typedef dynamic_matrix<scalar_type>                   matrix_type;

            static matrix_type impl(const mesh_type& msh, const cell_type& cl, const BQData& bqd)
            {
               const auto num_cell_dofs = howmany_dofs(bqd.cell_basis);

               matrix_type mass = matrix_type::Zero(num_cell_dofs, num_cell_dofs);

               const auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
               for (auto& qp : cell_quadpoints)
               {
                  const auto cphi = bqd.cell_basis.eval_functions(msh, cl, qp.point());
                  assert(cphi.size() == num_cell_dofs);

                  for(size_t j = 0; j < num_cell_dofs; j++){
                     for(size_t i = j; i < num_cell_dofs; i++){
                        mass(i,j) += qp.weight() * mm_prod(cphi[i], cphi[j]);
                     }
                  }
               }

               //lower part
               for(size_t j = 0; j < num_cell_dofs; j++){
                  for(size_t i = 0; i < j; i++){
                     mass(i,j) = mass(j,i);
                  }
               }

               return mass;
            }

            static matrix_type impl(const mesh_type& msh, const face_type& fc, const BQData& bqd)
            {
               const auto num_face_dofs = howmany_dofs(bqd.face_basis);

               matrix_type mass = matrix_type::Zero(num_face_dofs, num_face_dofs);

               const auto face_quadpoints = bqd.face_quadrature.integrate(msh, cl);
               for (auto& qp : face_quadpoints)
               {
                  const auto fphi = bqd.face_basis.eval_functions(msh, fc, qp.point());
                  assert(fphi.size() == num_face_dofs);

                  for(size_t j = 0; j < num_face_dofs; j++){
                     for(size_t i = j; i < num_face_dofs; i++){
                        mass(i,j) += qp.weight() * mm_prod(fphi[i], fphi[j]);
                     }
                  }
               }

               //lower part
               for(size_t j = 0; j < num_face_dofs; j++){
                  for(size_t i = 0; i < j; i++){
                     mass(i,j) = mass(j,i);
                  }
               }

               return mass;
            }
         };

      }  // namespace priv

      //cell mass matrix
      template<typename BQData>
      dynamic_matrix<typename BQData::mesh_type::scalar_type>
      mass_matrix(const typename BQData::mesh_type& msh, const typename BQData::mesh_type::cell_type& cl, const BQData& bqd)
      { return priv::mass_matrix_F<BQData, typename BQData::cell_basis_type>::impl(msh, cl, bqd); }

      //face mass matrix
      template<typename BQData>
      dynamic_matrix<typename BQData::mesh_type::scalar_type>
      mass_matrix(const typename BQData::mesh_type& msh, const typename BQData::mesh_type::face_type& fc, const BQData& bqd)
      { return priv::mass_matrix_F<BQData, typename BQData::cell_basis_type>::impl(msh, fc, bqd); }*/

   } // hho namespace
} // disk namespace