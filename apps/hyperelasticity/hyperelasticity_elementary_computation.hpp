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
#include "bases/bases_utils.hpp"
#include "mesh/mesh.hpp"
#include "bases/bases_ranges.hpp"
#include "hho/hho_nl_vector.hpp"
#include "BehaviorLaws/behaviorlaws.hpp"
#include "BehaviorLaws/maths_tensor.hpp"
#include "timecounter.h"

//#define USE_BLAS

namespace Hyperelasticity {


   template<typename Mesh, typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType,
   typename GradBasisType, typename GradQuadType>
   class Hyperelasticity
   {
      typedef Mesh                                mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;
      typedef typename mesh_type::face            face_type;

      typedef CellBasisType                       cell_basis_type;
      typedef CellQuadType                        cell_quadrature_type;
      typedef FaceBasisType                       face_basis_type;
      typedef FaceQuadType                        face_quadrature_type;
      typedef GradBasisType                       grad_basis_type;
      typedef GradQuadType                        grad_quadrature_type;

      typedef dynamic_matrix<scalar_type>         matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;

      cell_basis_type                             cell_basis;
      cell_quadrature_type                        cell_quadrature;

      face_basis_type                             face_basis;
      face_quadrature_type                        face_quadrature;

      grad_basis_type                             grad_basis;
      grad_quadrature_type                        grad_quadrature;

      size_t                                      m_degree;

      template<typename NeumannFunction>
      void
      add_NeumannConditions(const mesh_type& msh, const cell_type& cl, const NeumannFunction& g,
                            const std::vector<size_t>& boundary_neumann)
      {
         auto fcs = faces(msh, cl);

         size_t face_i = 0;
         for (auto& fc : fcs)
         {
            //Find if this face is a boundary face
            if(msh.is_boundary(fc)){
               //Find if this face is a boundary face with Neumann Condition
               if ( std::find(boundary_neumann.begin(), boundary_neumann.end(), msh.boundary_id(fc))
                   != boundary_neumann.end() ){
                  auto face_quadpoints = face_quadrature.integrate(msh, fc);
                  for (auto& qp : face_quadpoints)
                  {
                     auto fphi = face_basis.eval_functions(msh, fc, qp.point());

                     for(size_t i=0; i < fphi.size(); i++)
                        RTF(i) += qp.weight() * disk::mm_prod(g(qp.point()) , fphi[i]);
                  }
               }
            }
            face_i++;
         }

      }


   public:
      matrix_type     K_int;
      vector_type     RTF;
      double time_law, time_adapt_stab;
      scalar_type beta_adap;

      Hyperelasticity()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree);
         cell_quadrature     = cell_quadrature_type(2*m_degree);
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
         grad_basis          = grad_basis_type(m_degree);
         grad_quadrature     = grad_quadrature_type(2*m_degree);
      }

      Hyperelasticity(size_t degree)
      : m_degree(degree)
      {
         cell_basis          = cell_basis_type(m_degree);
         cell_quadrature     = cell_quadrature_type(2*(m_degree));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
         grad_basis          = grad_basis_type(m_degree);
         grad_quadrature     = grad_quadrature_type(2*m_degree);
      }

      matrix_type cell_mm;

      template<typename Function, typename NeumannFunction>
      void
      compute(const mesh_type& msh, const cell_type& cl, const Function& load, const NeumannFunction& neumann,
              const std::vector<size_t>& boundary_neumann, const matrix_type& GT,
              const vector_type& uTF, const ElasticityParameters elas_param,
              const bool adapt_stab = false)
      {
         time_law = 0.0;
         time_adapt_stab = 0.0;
         timecounter tc;
         const size_t DIM= msh.dimension;
         const size_t cpk = DIM * binomial(m_degree + DIM, m_degree);
         const size_t fpk = DIM * binomial(m_degree  + DIM -1, m_degree);
         const size_t gpk = DIM * DIM * binomial(m_degree + DIM, m_degree);

         auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();

         const size_t num_grad_dofs = DIM * grad_basis.range(0, m_degree).size();
         const size_t num_cell_dofs = cell_basis.range(0, m_degree).size();
         const size_t num_face_dofs = face_basis.size();
         const size_t num_total_dofs = num_cell_dofs + num_faces * num_face_dofs;

         assert(num_grad_dofs == gpk);
         assert(num_cell_dofs == cpk);
         assert(num_face_dofs == fpk);


         matrix_type AT = matrix_type::Zero(num_grad_dofs, num_grad_dofs);
         vector_type aT = vector_type::Zero(num_grad_dofs);

         RTF = vector_type::Zero(num_total_dofs);

         vector_type GT_uTF = GT * uTF;

         NeoHookeanLaw<scalar_type>  law(elas_param.mu, elas_param.lambda, elas_param.type_law);

         auto grad_quadpoints = grad_quadrature.integrate(msh, cl);

         beta_adap = elas_param.tau;

         for (auto& qp : grad_quadpoints)
         {
            auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
            auto gphi = grad_basis.eval_functions(msh, cl, qp.point());
            assert(num_grad_dofs == gphi.size());
            assert(num_cell_dofs == c_phi.size());

            // Compute local gradient and norm
            auto GT_iqn = disk::compute_gradient_matrix_pt(GT_uTF, gphi);
            auto FT_iqn = compute_FTensor(GT_iqn);

            //Compute bahavior
            tc.tic();
            auto tensor_behavior = law.compute_whole_PK1(FT_iqn);
            tc.toc();
            time_law += tc.to_double();

            if(adapt_stab)
            {
               tc.tic();
               Eigen::SelfAdjointEigenSolver<decltype(tensor_behavior.second)> es;
               decltype(tensor_behavior.second) Arow = changeFormatRowTensor(tensor_behavior.second);
               es.compute(Arow);
               scalar_type ev_min = es.eigenvalues().minCoeff();
               //std::cout << ev_min << std::endl;
               beta_adap = std::max(elas_param.tau, std::max( beta_adap, std::abs(ev_min)));
               time_adapt_stab += tc.to_double();
            }


            for(std::size_t i = 0; i < num_grad_dofs; i++) {
               auto Agphi_i = tm_prod(tensor_behavior.second, gphi[i]);
               for(std::size_t j = i; j < num_grad_dofs; j++) {
                  //compute (Gkt v, A(u) : Gkt du)
                  AT(i,j) += qp.weight() * disk::mm_prod(Agphi_i, gphi[j]);
               }
               // compute (PK1(u), G^k_T v)_T
               aT(i) += qp.weight() * disk::mm_prod(tensor_behavior.first, gphi[i]);
            }

            //compute (f,v)_T
            for(std::size_t i = 0; i < num_cell_dofs; i++) {
               RTF(i) += qp.weight() * disk::mm_prod(load(qp.point()) , c_phi[i]);
            }
         }

         //lower part AT
         for(std::size_t i = 1; i < num_grad_dofs; i++)
            for(std::size_t j = 0; j < i; j++)
               AT(i,j) = AT(j,i);

         K_int = GT.transpose() * AT * GT;
         RTF -= GT.transpose() * aT;

         //Add Neumann condition
         add_NeumannConditions(msh, cl, neumann, boundary_neumann);

         assert(K_int.rows() == num_total_dofs);
         assert(K_int.cols() == num_total_dofs);
         assert(RTF.rows() == num_total_dofs);
      }

   template<typename Function, typename NeumannFunction>
   void
   compute2(const mesh_type& msh, const cell_type& cl, const Function& load, const NeumannFunction& neumann,
           const std::vector<size_t>& boundary_neumann, const matrix_type& GT,
           const vector_type& uTF, const ElasticityParameters elas_param)
   {
      time_law = 0.0;
      time_adapt_stab = 0.0;
      timecounter tc;
      const size_t DIM= msh.dimension;
      const size_t cpk = DIM * binomial(m_degree + DIM, m_degree);
      const size_t fpk = DIM * binomial(m_degree  + DIM -1, m_degree);
      const size_t gpk = DIM * binomial(m_degree + 1 + DIM, m_degree + 1);

      auto fcs = faces(msh, cl);
      const size_t num_faces = fcs.size();

      auto grad_range = cell_basis.range(1, m_degree+1);

      const size_t num_grad_dofs = grad_range.size();
      const size_t num_cell_dofs = cell_basis.range(0, m_degree).size();
      const size_t num_face_dofs = face_basis.size();
      const size_t num_total_dofs = num_cell_dofs + num_faces * num_face_dofs;

      assert(num_grad_dofs == gpk);
      assert(num_cell_dofs == cpk);
      assert(num_face_dofs == fpk);

      matrix_type AT = matrix_type::Zero(num_grad_dofs, num_grad_dofs);
      vector_type aT = vector_type::Zero(num_grad_dofs);

      RTF = vector_type::Zero(num_total_dofs);

      vector_type GT_uTF = GT * uTF;

      NeoHookeanLaw<scalar_type>  law(elas_param.mu, elas_param.lambda, elas_param.type_law);

      auto cell_quadpoints = cell_quadrature.integrate(msh, cl);

      scalar_type beta_adap = elas_param.tau;

      for (auto& qp : cell_quadpoints)
      {
         auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
         auto g_phi = cell_basis.eval_gradients(msh, cl, qp.point());
         // delete DIM firs element which are zero
         g_phi.erase(g_phi.begin(), g_phi.begin() + grad_range.from());
         assert(num_grad_dofs == g_phi.size());
         assert(num_cell_dofs == c_phi.size());

         // Compute local gradient and norm
         auto GT_iqn = disk::compute_gradient_matrix_pt(GT_uTF, g_phi);
         auto FT_iqn = compute_FTensor(GT_iqn);

         //Compute bahavior
         tc.tic();
         auto tensor_behavior = law.compute_whole_PK1(FT_iqn);
         tc.toc();
         time_law += tc.to_double();

         if(elas_param.adaptative_stab)
         {
            tc.tic();
            Eigen::SelfAdjointEigenSolver<decltype(tensor_behavior.second)> es;
            decltype(tensor_behavior.second) Arow = changeFormatRowTensor(tensor_behavior.second);
            es.compute(Arow);
            scalar_type ev_min = es.eigenvalues().minCoeff();
            std::cout << ev_min << std::endl;
            beta_adap = std::max(elas_param.tau, std::max( beta_adap, std::abs(ev_min)));
            time_adapt_stab += tc.to_double();
         }


         for(std::size_t i = 0; i < num_grad_dofs; i++) {
            auto Agphi_i = tm_prod(tensor_behavior.second, g_phi[i]);
            for(std::size_t j = i; j < num_grad_dofs; j++) {
               //compute (Gkt v, A(u) : Gkt du)
               AT(i,j) += qp.weight() * disk::mm_prod(Agphi_i, g_phi[j]);
            }
            // compute (PK1(u), G^k_T v)_T
            aT(i) += qp.weight() * disk::mm_prod(tensor_behavior.first, g_phi[i]);
         }

         //compute (f,v)_T
         for(std::size_t i = 0; i < c_phi.size(); i++) {
            RTF(i) += qp.weight() * disk::mm_prod(load(qp.point()) , c_phi[i]);
         }
      }

      //lower part AT
      for(std::size_t i = 1; i < num_grad_dofs; i++)
         for(std::size_t j = 0; j < i; j++)
            AT(i,j) = AT(j,i);

      K_int = GT.transpose() * AT * GT;
      RTF -= GT.transpose() * aT;

      //Add Neumann condition
      add_NeumannConditions(msh, cl, neumann, boundary_neumann);

      assert(K_int.rows() == num_total_dofs);
      assert(K_int.cols() == num_total_dofs);
      assert(RTF.rows() == num_total_dofs);
   }

};

}//end namespace hyperelasticity
