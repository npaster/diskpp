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
#include "bases/bases_utils.hpp"
#include "mesh/mesh.hpp"
#include "hho/hho_utils.hpp"
#include "BehaviorLaws/behaviorlaws.hpp"
#include "BehaviorLaws/maths_tensor.hpp"
#include "mechanics/deformation_tensors.hpp"
#include "timecounter.h"

//#define USE_BLAS

namespace Hyperelasticity {

//    namespace priv {
//       template< typename GradBasis >
//       struct poly_size_F {
//          static size_t impl(const size_t DIM, const size_t grad_degree)
//          { throw std::invalid_argument("Unknown GradBasis"); }
//       };
//
//       template<>
//       struct poly_size_F<disk::Raviart_Thomas_matrix_basis> {
//          static size_t impl(const size_t DIM, const size_t grad_degree)
//          { return DIM * DIM * binomial(grad_degree-1 + DIM, grad_degree-1); }
//       };
//
//       template<>
//       struct poly_size_F<disk::scaled_monomial_matrix_basis> {
//          static size_t impl(const size_t DIM, const size_t grad_degree)
//          { return DIM * DIM * binomial(grad_degree + DIM, grad_degree); }
//       };
//
//   }  // namespace priv


   template<typename BQData>
   class Hyperelasticity
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;

      typedef dynamic_matrix<scalar_type>         matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;

      const BQData&                               m_bqd;


      // Optimzation of the product A:G
      //attention grad_degree = 0

//       template<typename GradBasis>
//       size_t poly_size(const size_t DIM, const size_t grad_degree)
//       { return priv::poly_size_F<GradBasis>::impl(const size_t DIM, const size_t grad_degree); }

      template<int DIM>
      std::vector<static_matrix<scalar_type, DIM, DIM>>
      compute_A_gphi(const static_tensor<scalar_type, DIM>& tens,
                     const std::vector<static_matrix<scalar_type, DIM, DIM>>& gphi ) const
      {
         const size_t grad_basis_size = gphi.size();
         const size_t grad_degree = m_bqd.grad_degree();
         const size_t poly_space = DIM * DIM * binomial(grad_degree-1 + DIM, grad_degree-1);
         const size_t DIM2 = DIM * DIM;

         std::vector<static_matrix<scalar_type, DIM, DIM>> Aphi;
         Aphi.reserve(grad_basis_size);

         //poly classique
         for(std::size_t i = 0; i < poly_space; i += DIM2) {
            size_t row = i;
            for(size_t k = 0; k < DIM; k++ ){//depend de l'ordre des bases
               for(size_t l = 0; l < DIM; l++ ){//depend de l'ordre des bases
                  Aphi.push_back(tm_prod(tens, gphi[row], l, k));
                  row++;
               }
            }
         }

         // RT space
         for(std::size_t i = poly_space; i < grad_basis_size; i ++) {
            Aphi.push_back(tm_prod(tens, gphi[i]));
         }

         return Aphi;
      }



      template<int DIM>
      void
      compute_A_gphi_gphi(const static_tensor<scalar_type, DIM>& A,
                          const std::vector<static_matrix<scalar_type, DIM, DIM>>& gphi,
                          const scalar_type weight, matrix_type& AT) const
      {
         const size_t grad_basis_size = gphi.size();
         const size_t grad_degree = m_bqd.grad_degree();
         const size_t poly_space = DIM * DIM * binomial(grad_degree-1 + DIM, grad_degree-1);
         const size_t DIM2 = DIM * DIM;

         assert(AT.rows() == grad_basis_size);
         assert(AT.cols() == grad_basis_size);

         // compute A:gphi

         const auto A_gphi = compute_A_gphi(A, gphi);
         assert(grad_basis_size == A_gphi.size());

         //espace poly classique

         for(size_t j = 0; j < poly_space; j += DIM2 ){
            size_t col = j;
            for(size_t k = 0; k < DIM; k++ ){//depend de l'ordre des bases
               for(size_t l = 0; l < DIM; l++ ){//depend de l'ordre des bases
                  for(size_t i = col; i < poly_space; i++){
                     AT(i,col) += weight * A_gphi[i](l,k) * gphi[col](l,k);
                  }
                  col++;
               }
            }
         }

         //espace RT

         for(std::size_t i = poly_space; i < grad_basis_size; i ++) {
            for(std::size_t j = 0; j < grad_basis_size; j ++) {
               AT(i,j) += weight * disk::mm_prod(A_gphi[i], gphi[j]);
            }
         }
      }


      template<int DIM>
      void
      compute_PK1_gphi(const static_matrix<scalar_type, DIM, DIM>& PK1,
                       const std::vector<static_matrix<scalar_type, DIM, DIM>>& gphi,
                       const scalar_type weight, vector_type& aT) const
      {
         const size_t grad_basis_size = gphi.size();
         const size_t grad_degree = m_bqd.grad_degree();
         const size_t poly_space = DIM * DIM * binomial(grad_degree-1 + DIM, grad_degree-1);
         const size_t DIM2 = DIM * DIM;

         assert(aT.size() == grad_basis_size);

         //espace poly classique

         for(std::size_t i = 0; i < poly_space; i+= DIM2) {
            size_t row = i;
            for(size_t k = 0; k < DIM; k++ ){//depend de l'ordre des bases
               for(size_t l = 0; l < DIM; l++ ){//depend de l'ordre des bases
                  // compute (PK1(u), G^k_T v)_T
                  aT(row) += weight * PK1(l,k) * gphi[row](l,k);
                  row++;
               }
            }
         }

         //espace RT

         for(std::size_t i = poly_space; i < grad_basis_size; i++) {
            aT(i) += weight * disk::mm_prod(PK1, gphi[i]);
         }
      }



   public:
      matrix_type     K_int;
      vector_type     RTF;
      double time_law, time_adapt_stab;

      Hyperelasticity(const BQData& bqd) : m_bqd(bqd)
      {}

      template<typename Function>
      void
      compute(const mesh_type& msh, const cell_type& cl, const Function& load,
              const matrix_type& GT, const vector_type& uTF, const ElasticityParameters elas_param)
      {
         const size_t DIM= msh.dimension;
         const size_t cell_degree = m_bqd.cell_degree();
         const size_t face_degree = m_bqd.face_degree();
         const size_t grad_degree = m_bqd.grad_degree();
         const size_t cell_basis_size = (m_bqd.cell_basis.range(0, cell_degree)).size();
         const size_t face_basis_size = (m_bqd.face_basis.range(0, face_degree)).size();
         const size_t grad_basis_size = (m_bqd.grad_basis.range(0, grad_degree)).size();

         time_law = 0.0;
         time_adapt_stab = 0.0;
         timecounter tc;

         const auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();
         const size_t num_total_dofs = cell_basis_size + num_faces * face_basis_size;

         matrix_type AT = matrix_type::Zero(grad_basis_size, grad_basis_size);
         vector_type aT = vector_type::Zero(grad_basis_size);

         RTF = vector_type::Zero(num_total_dofs);

         const vector_type GT_uTF = GT * uTF;

         const NeoHookeanLaw<scalar_type>  law(elas_param.mu, elas_param.lambda, elas_param.type_law);
         //CavitationLaw<scalar_type>  law(elas_param.mu, elas_param.lambda, elas_param.type_law);

         const auto grad_quadpoints = m_bqd.grad_quadrature.integrate(msh, cl);

         for (auto& qp : grad_quadpoints)
         {
            const auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());
            assert(grad_basis_size == gphi.size());

            // Compute local gradient and norm
            const auto GT_iqn = disk::hho::eval_gradient(GT_uTF, gphi);
            const auto FT_iqn = disk::mechanics::convertGtoF(GT_iqn);

            //Compute bahavior
            tc.tic();
            const auto tensor_behavior = law.compute_whole_PK1(FT_iqn);
            tc.toc();
            time_law += tc.to_double();

            // compute (A(u):G^k_T du, G^k_T v)_T
            this->compute_A_gphi_gphi(tensor_behavior.second, gphi, qp.weight(), AT);

            // compute (PK1(u), G^k_T v)_T
            this->compute_PK1_gphi(tensor_behavior.first, gphi, qp.weight(), aT);
         }


         const auto cell_quadpoints = m_bqd.cell_quadrature.integrate(msh, cl);

         for (auto& qp : cell_quadpoints)
         {
            const auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point());
            //compute (f,v)_T
            const auto load_qp = load(qp.point());
            for(std::size_t i = 0; i < cell_basis_size; i += DIM) {
               size_t row = i;
               for(size_t k = 0; k < DIM; k++ ){//depend de l'ordre des bases
                  RTF(row) += qp.weight() * load_qp(k) * c_phi[row](k);
                  row++;
               }
            }
         }

         //lower part AT
         for(std::size_t i = 0; i < grad_basis_size; i++)
            for(std::size_t j = i; j < grad_basis_size; j++)
               AT(i,j) = AT(j,i);

         K_int = GT.transpose() * AT * GT;
         RTF -= GT.transpose() * aT;

         assert(K_int.rows() == num_total_dofs);
         assert(K_int.cols() == num_total_dofs);
         assert(RTF.rows() == num_total_dofs);
      }
};

}//end namespace hyperelasticity
