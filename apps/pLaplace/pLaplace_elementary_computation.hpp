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
#include "bases/bases_ranges.hpp"
#include "BehaviorLaws/behaviorlaws.hpp"
#include "BehaviorLaws/maths_tensor.hpp"
#include "timecounter.h"

//#define USE_BLAS
#define FILL_COLMAJOR

namespace pLaplace {
   
   
   //version provisoire à compléter
   template<typename Mesh, typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType>
   class plaplace_like_stabilization
   {
      typedef Mesh                            mesh_type;
      typedef typename mesh_type::scalar_type scalar_type;
      typedef typename mesh_type::cell        cell_type;
      typedef typename mesh_type::face        face_type;
      
      typedef CellBasisType                   cell_basis_type;
      typedef CellQuadType                    cell_quadrature_type;
      typedef FaceBasisType                   face_basis_type;
      typedef FaceQuadType                    face_quadrature_type;
      
      typedef dynamic_matrix<scalar_type>     matrix_type;
      typedef dynamic_vector<scalar_type>     vector_type;
      
      cell_basis_type                         cell_basis;
      cell_quadrature_type                    cell_quadrature;
      
      face_basis_type                         face_basis;
      face_quadrature_type                    face_quadrature;
      
      size_t                                  m_degree;
      
   public:
      matrix_type     data;
      
      plaplace_like_stabilization()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }
      
      plaplace_like_stabilization(size_t degree)
      : m_degree(degree)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }
      
      void compute(const mesh_type& msh, const cell_type& cl, const matrix_type& gradrec_oper)
      {
         matrix_type mass_mat = matrix_type::Zero(cell_basis.size(), cell_basis.size());
         
         auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
            mass_mat += qp.weight() * c_phi * c_phi.transpose();
         }
         
         auto zero_range         = cell_basis.range(0, m_degree);
         auto one_range          = cell_basis.range(1, m_degree+1);
         
         // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)
         
         //Step 1: compute \pi_T^k p_T^k v (third term).
         matrix_type M1 = take(mass_mat, zero_range, zero_range);
         matrix_type M2 = take(mass_mat, zero_range, one_range);
         matrix_type proj1 = -M1.llt().solve(M2*gradrec_oper);
         
         //Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
         matrix_type I_T = matrix_type::Identity(zero_range.size(), zero_range.size());
         proj1.block(0, 0, zero_range.size(), zero_range.size()) += I_T;
         
         auto fcs = faces(msh, cl);
         size_t num_faces = fcs.size();
         
         size_t num_cell_dofs = cell_basis.range(0, m_degree).size();
         size_t num_face_dofs = face_basis.size();
         
         disk::dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);
         
         data = matrix_type::Zero(dsr.total_size(), dsr.total_size());
         
         // Step 3: project on faces (eqn. 21)
         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto current_face_range = dsr.face_range(face_i);
            size_t fbs = current_face_range.size();
            auto h = diameter(msh, /*fcs[face_i]*/cl);
            auto fc = fcs[face_i];
            
            matrix_type face_mass_matrix    = matrix_type::Zero(fbs, fbs);
            matrix_type face_trace_matrix   = matrix_type::Zero(fbs, cell_basis.size());
            
            auto face_quadpoints = face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
               auto f_phi = face_basis.eval_functions(msh, fc, qp.point());
               auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
               auto q_f_phi = qp.weight() * f_phi;
               face_mass_matrix += q_f_phi * f_phi.transpose();
               face_trace_matrix += q_f_phi * c_phi.transpose();
            }
            
            Eigen::LLT<matrix_type> piKF;
            piKF.compute(face_mass_matrix);
            
            // Step 3a: \pi_F^k( v_F - p_T^k v )
            auto face_range = current_face_range.remove_offset();
            matrix_type MR1 = take(face_trace_matrix, face_range, one_range);
            matrix_type proj2 = piKF.solve(MR1*gradrec_oper);
            matrix_type I_F = matrix_type::Identity(fbs, fbs);
            auto block_offset = current_face_range.min();
            proj2.block(0, block_offset, fbs, fbs) -= I_F;
            
            // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
            matrix_type MR2 = take(face_trace_matrix, face_range, zero_range);
            matrix_type proj3 = piKF.solve(MR2*proj1);
            
            matrix_type BRF = proj2 + proj3;
            
            data += BRF.transpose() * face_mass_matrix * BRF / h;
         }
      }
   };
   
   template<typename Mesh, typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType,
   typename GradBasisType, typename GradQuadType>
   class pLaplacian
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
      
   public:
      matrix_type     K_int;
      vector_type     RTF;      

      
      pLaplacian()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree);
         cell_quadrature     = cell_quadrature_type(2*m_degree);
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
         grad_basis          = grad_basis_type(m_degree);
         grad_quadrature     = grad_quadrature_type(2*m_degree);
      }
      
      pLaplacian(size_t degree)
      : m_degree(degree)
      {
         cell_basis          = cell_basis_type(m_degree);
         cell_quadrature     = cell_quadrature_type(2*m_degree);
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
         grad_basis          = grad_basis_type(m_degree);
         grad_quadrature     = grad_quadrature_type(2*m_degree);
      }
      
      matrix_type cell_mm;
      
      template<typename Function>
      void
      compute(const mesh_type& msh, const cell_type& cl, const Function& load, const matrix_type& GT,
               const vector_type& uTF, const size_t p)
      {
         auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();
         
         const size_t num_grad_dofs = grad_basis.range(0, m_degree).size();
         const size_t num_cell_dofs = cell_basis.range(0, m_degree).size();
         const size_t num_face_dofs = face_basis.size();
         const size_t num_total_dofs = num_cell_dofs + num_faces * num_face_dofs;
         
         
         matrix_type AT = matrix_type::Zero(num_grad_dofs, num_grad_dofs);
         vector_type aT = vector_type::Zero(num_grad_dofs);
         
         RTF = vector_type::Zero(num_total_dofs);
         
         vector_type GT_uTF = GT * uTF;
         
         pLaplaceLaw<scalar_type>  law(p);
         
         auto grad_quadpoints = grad_quadrature.integrate(msh, cl);
         
         for (auto& qp : grad_quadpoints)
         {
            matrix_type c_phi = cell_basis.eval_functions(msh, cl, qp.point(), 0, m_degree);
            auto gphi = grad_basis.eval_functions(msh, cl, qp.point());
            assert(num_grad_dofs == gphi.size());
            
            // Compute local gradient and norm
            auto GT_iqn = disk::compute_gradient_vector_pt<grad_basis_type>(msh, cl, GT_uTF, qp.point(), m_degree);
            auto tensor_behavior = law.compute_whole(GT_iqn);
            
            for(std::size_t i = 0; i < num_grad_dofs; i++) {
               for(std::size_t j = 0; j < num_grad_dofs; j++) {
                  AT(i,j) += qp.weight() * disk::mm_prod(gphi.at(i), disk::mm_prod(tensor_behavior.second, gphi.at(j)));
               } // for j
               aT(i) += qp.weight() * disk::mm_prod(tensor_behavior.first, gphi.at(i));
            } // for i
            
            for(std::size_t i = 0; i < num_cell_dofs; i++) {
               RTF(i) += qp.weight() * load(qp.point()) * c_phi(i);
            } // for i
         }
         
         K_int = GT.transpose() * AT * GT;
         RTF -= GT.transpose() * aT;
         
         assert(K_int.rows() == num_total_dofs);
         assert(K_int.cols() == num_total_dofs);
         assert(RTF.rows() == num_total_dofs);
      }

   };
   
}