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
#include "hho/hho_nl.hpp"
#include "BehaviorLaws/behaviorlaws.hpp"
#include "timecounter.h"

//#define USE_BLAS
#define FILL_COLMAJOR

namespace disk {
   
   template<typename GradBasisType, typename Mesh, size_t DIM>
   static_matrix<typename Mesh::scalar_type,DIM, DIM>
   compute_gradient_matrix_pt(const Mesh& msh, const typename Mesh::cell& cl,
                              const dynamic_vector<typename Mesh::scalar_type>& gradrec_coeff,
                              const point<typename Mesh::scalar_type,DIM>& pt, const size_t degree)
   {
      typedef typename Mesh::scalar_type               scalar_type;
      typedef static_matrix<scalar_type,DIM, DIM>      gradient_value_type;
      
      GradBasisType grad_basis          = GradBasisType(degree);
      
      gradient_value_type ret = gradient_value_type::Zero();
      
      auto gphi = grad_basis.eval_functions(msh, cl, pt);
      
      for(size_t i = 0; i < gphi.size(); i++){
         ret += gradrec_coeff(i) * gphi.at(i);
      }
      
      return ret;
   }
   
   
   template<typename Mesh,typename FaceBasisType, typename FaceQuadType>
   class assembler_nl_vector
   {
      typedef Mesh                                mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;
      typedef typename mesh_type::face            face_type;
      
      typedef FaceBasisType                       face_basis_type;
      typedef FaceQuadType                        face_quadrature_type;
      
      typedef dynamic_matrix<scalar_type>         matrix_type;
      
      face_basis_type                             face_basis;
      face_quadrature_type                        face_quadrature;
      
      size_t                                      m_degree;
      
      typedef Eigen::Triplet<scalar_type>         triplet_type;
      
      std::vector<triplet_type>                   m_triplets;
      size_t                                      m_num_unknowns;
      
   public:
      
      typedef Eigen::SparseMatrix<scalar_type>    sparse_matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;
      
      sparse_matrix_type      matrix;
      vector_type             rhs;
      
      assembler_nl_vector()                 = delete;
      
      assembler_nl_vector(const mesh_type& msh, size_t degree)
      : m_degree(degree)
      {
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
         
         m_num_unknowns = face_basis.size() * (msh.faces_size() + msh.boundary_faces_size());
         matrix = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
         rhs = vector_type::Zero(m_num_unknowns);
      }
      
      template<typename LocalContrib>
      void
      assemble(const mesh_type& msh, const cell_type& cl, const LocalContrib& lc)
      {
         auto fcs = faces(msh, cl);
         std::vector<size_t> l2g(fcs.size() * face_basis.size());
         for (size_t face_i = 0; face_i < fcs.size(); face_i++)
         {
            auto fc = fcs[face_i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
               throw std::invalid_argument("This is a bug: face not found");
            
            auto face_id = eid.second;
            
            auto face_offset = face_id * face_basis.size();
            
            auto pos = face_i * face_basis.size();
            
            for (size_t i = 0; i < face_basis.size(); i++)
               l2g[pos+i] = face_offset+i;
         }
         
         assert(lc.first.rows() == lc.first.cols());
         assert(lc.first.rows() == lc.second.size());
         assert(lc.second.size() == l2g.size());
         
         //std::cout << lc.second.size() << " " << l2g.size() << std::endl;
         
         #ifdef FILL_COLMAJOR
         for (size_t j = 0; j < lc.first.cols(); j++)
         {
            for (size_t i = 0; i < lc.first.rows(); i++)
               m_triplets.push_back( triplet_type( l2g.at(i), l2g.at(j), lc.first(i,j) ) );
            
            rhs(l2g.at(j)) += lc.second(j);
         }
         #else
         for (size_t i = 0; i < lc.first.rows(); i++)
         {
            for (size_t j = 0; j < lc.first.cols(); j++)
               m_triplets.push_back( triplet_type( l2g.at(i), l2g.at(j), lc.first(i,j) ) );
            
            rhs(l2g.at(i)) += lc.second(i);
         }
         #endif
      }
      
      template<typename Function>
      void
      impose_boundary_conditions(const mesh_type& msh, const Function& bc, const std::vector<vector_type>& sol_faces,
                                 const std::vector<vector_type>& sol_lagr)
      {
         size_t fbs = face_basis.size();
         size_t face_i = 0;
         for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
         {
            auto bfc = *itor;
            
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first)
               throw std::invalid_argument("This is a bug: face not found");
            
            auto face_id = eid.second;
            
            auto face_offset = face_id * fbs;
            auto face_offset_lagrange = (msh.faces_size() + face_i) * fbs;
            
            auto fqd = face_quadrature.integrate(msh, bfc);
            
            matrix_type MFF     = matrix_type::Zero(fbs, fbs);
            vector_type rhs_f   = vector_type::Zero(fbs);
            
            for (auto& qp : fqd)
            {
               auto f_phi = face_basis.eval_functions(msh, bfc, qp.point());
               for(size_t i = 0; i < fbs; i++)
                  for(size_t j = 0; j < fbs; j++)
                     MFF(i,j) += qp.weight() * mm_prod(f_phi.at(i), f_phi.at(j));
                  
               for(size_t i = 0; i < fbs; i++)
                  rhs_f(i) += qp.weight() * mm_prod( f_phi.at(i),  bc(qp.point()));
            }
            
            
            rhs_f -= MFF * sol_faces.at(face_id);
            
            vector_type rhs_l = MFF * sol_lagr.at(face_i);
            
            #ifdef FILL_COLMAJOR
            for (size_t j = 0; j < MFF.cols(); j++)
            {
               for (size_t i = 0; i < MFF.rows(); i++)
               {
                  m_triplets.push_back( triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i,j)) );
                  m_triplets.push_back( triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i,j)) );
               }
               rhs(face_offset_lagrange+j) = rhs_f(j);
               rhs(face_offset+j) -= rhs_l(j);
            }
            #else
            for (size_t i = 0; i < MFF.rows(); i++)
            {
               for (size_t j = 0; j < MFF.cols(); j++)
               {
                  m_triplets.push_back( triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i,j)) );
                  m_triplets.push_back( triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i,j)) );
               }
               rhs(face_offset_lagrange+i) = rhs_f(i);
               rhs(face_offset+j) -= rhs_l(j);
            }
            #endif
            
            face_i++;
         }
      }
      
      void
      finalize()
      {
         matrix.setFromTriplets(m_triplets.begin(), m_triplets.end());
         m_triplets.clear();
      }
      
      void
      finalize(sparse_matrix_type& mat, vector_type& vec)
      {
         mat = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
         mat.setFromTriplets(m_triplets.begin(), m_triplets.end());
         m_triplets.clear();
         vec = rhs;
      }
   };
   
   
   
   
   template<typename Mesh, typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType,
   typename GradBasisType, typename GradQuadType>
   class Laplacian
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

      
      Laplacian()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree);
         cell_quadrature     = cell_quadrature_type(2*m_degree);
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
         grad_basis          = grad_basis_type(m_degree);
         grad_quadrature     = grad_quadrature_type(2*m_degree);
      }
      
      Laplacian(size_t degree)
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
         
         LinearElasticityLaw<scalar_type>  law(10.0);
         
         auto grad_quadpoints = grad_quadrature.integrate(msh, cl);
         
         for (auto& qp : grad_quadpoints)
         {
            auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
            auto gphi = grad_basis.eval_functions(msh, cl, qp.point());
            assert(num_grad_dofs == gphi.size());
            assert(num_cell_dofs == c_phi.size());
            
            // Compute local gradient and norm
            auto GT_iqn = compute_gradient_matrix_pt<grad_basis_type>(msh, cl, GT_uTF, qp.point(), m_degree);
            auto tensor_behavior = law.compute_whole_PK1(GT_iqn);
            
            for(std::size_t i = 0; i < num_grad_dofs; i++) {
               for(std::size_t j = 0; j < num_grad_dofs; j++) {
                  AT(i,j) += qp.weight() * mm_prod(gphi.at(i), tm_prod(tensor_behavior.second, gphi.at(j)));
               } // for j
               aT(i) += qp.weight() * mm_prod(tensor_behavior.first, gphi.at(i));
            } // for i
            
            for(std::size_t i = 0; i < num_cell_dofs; i++) {
               RTF(i) += qp.weight() * mm_prod(load(qp.point()) , c_phi.at(i));
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