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
#include "timecounter.h"

//#define USE_BLAS
#define FILL_COLMAJOR

namespace disk {
   
   template<typename GradBasisType, typename Mesh, size_t DIM>
   static_vector<typename Mesh::scalar_type,DIM>
   compute_gradient_vector_pt(const Mesh& msh, const typename Mesh::cell& cl,
                              const dynamic_vector<typename Mesh::scalar_type>& gradrec_coeff,
                              const point<typename Mesh::scalar_type,DIM>& pt, const size_t degree)
   {
      typedef typename Mesh::scalar_type               scalar_type;
      typedef static_vector<scalar_type,DIM>            gradient_value_type;
      
      GradBasisType grad_basis          = GradBasisType(degree);
      
      gradient_value_type ret = gradient_value_type::Zero(DIM);
      
      auto gphi = grad_basis.eval_functions(msh, cl, pt);
      
      for(size_t i = 0; i < gphi.size(); i++){
         ret += gradrec_coeff(i) * gphi.at(i);
      }
      
      return ret;
   }
   
   
   template<typename Mesh,typename FaceBasisType, typename FaceQuadType>
   class assembler_nl_scalar
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
      
      assembler_nl_scalar()                 = delete;
      
      assembler_nl_scalar(const mesh_type& msh, size_t degree)
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
               MFF += qp.weight() * f_phi * f_phi.transpose();
               rhs_f += qp.weight() * f_phi * bc(qp.point());
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
}