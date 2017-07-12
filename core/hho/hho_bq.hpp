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
#include "timecounter.h"
//#include "contrib/sol2/sol.hpp"

//#define USE_BLAS
#define FILL_COLMAJOR

namespace disk {

   template<typename Mesh,
   template<typename, typename> class BasisFunction,
   template<typename, typename> class BasisGradient,
   template<typename, typename> class Quadrature>
   class basis_quadrature_data_full /* this name really sucks */
   {
   public:
      typedef Mesh                            mesh_type;
      typedef typename mesh_type::cell        cell_type;
      typedef typename mesh_type::face        face_type;

      typedef BasisFunction<mesh_type, cell_type>         cell_basis_type;
      typedef BasisFunction<mesh_type, face_type>         face_basis_type;
      typedef BasisGradient<mesh_type, cell_type>         grad_basis_type;

      typedef Quadrature<mesh_type, cell_type>    cell_quad_type;
      typedef Quadrature<mesh_type, face_type>    face_quad_type;
      typedef Quadrature<mesh_type, cell_type>    grad_quad_type;

      cell_basis_type     cell_basis;
      face_basis_type     face_basis;
      cell_quad_type      cell_quadrature;
      face_quad_type      face_quadrature;
      grad_basis_type     grad_basis;
      grad_quad_type      grad_quadrature;

   private:
      size_t  m_cell_degree, m_face_degree, m_grad_degree;

      void init(void)
      {
         cell_basis          = cell_basis_type(m_cell_degree);
         face_basis          = face_basis_type(m_face_degree);
         grad_basis          = grad_basis_type(m_grad_degree);
         cell_quadrature     = cell_quad_type(2 * m_cell_degree);
         face_quadrature     = face_quad_type(m_face_degree + m_grad_degree);
         grad_quadrature     = grad_quad_type(2 * m_grad_degree);
      }

   public:
      basis_quadrature_data_full() : m_cell_degree(1), m_face_degree(1), m_grad_degree(1)
      {
         init();
      }

      basis_quadrature_data_full(const size_t cell_degree, const size_t face_degree, const size_t grad_degree)
      {
         if ( (cell_degree + 1 < face_degree) or (cell_degree > face_degree + 1) )
            throw std::invalid_argument("Invalid cell degree");

         if ( (cell_degree + 1 < grad_degree) or (grad_degree < cell_degree) )
            throw std::invalid_argument("Invalid grad degree");

         m_cell_degree = cell_degree;
         m_face_degree = face_degree;
         m_grad_degree = grad_degree;

         init();
      }

      size_t cell_degree(void) const { return m_cell_degree; }
      size_t face_degree(void) const { return m_face_degree; }
      size_t grad_degree(void) const { return m_grad_degree; }
   };


   template<typename BQData>
   class gradient_reconstruction_full_bq
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;

      typedef dynamic_matrix<scalar_type>         matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;

      const BQData&                               m_bqd;

      matrix_type     m_oper;
      matrix_type     m_data;

   public:

      matrix_type oper() const {return m_oper;}
      matrix_type data() const {return m_data;}

      gradient_reconstruction_full_bq(const BQData& bqd) : m_bqd(bqd)
      {}

      void compute(const mesh_type& msh, const cell_type& cl, const bool compute_data = true)
      {
         const size_t DIM = msh.dimension;
         const size_t cell_degree = m_bqd.cell_degree();
         const size_t face_degree = m_bqd.face_degree();
         const size_t grad_degree = m_bqd.grad_degree();
         const size_t num_cell_dofs = (m_bqd.cell_basis.range(0, cell_degree)).size();
         const size_t num_face_dofs = (m_bqd.face_basis.range(0, face_degree)).size();
         const size_t num_grad_dofs = (m_bqd.grad_basis.range(0, grad_degree)).size();


         auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();

         assert(num_grad_dofs == DIM * binomial( grad_degree +DIM, grad_degree));

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         assert(dsr.total_size() == (num_cell_dofs + num_faces * num_face_dofs));

         matrix_type BG = matrix_type::Zero(num_grad_dofs, dsr.total_size());

         matrix_type MG = matrix_type::Zero(num_grad_dofs, num_grad_dofs);


         auto grad_quadpoints = m_bqd.grad_quadrature.integrate(msh, cl);
         for (auto& qp : grad_quadpoints)
         {
            auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());
            assert(num_grad_dofs == gphi.size());

            auto dphi = m_bqd.cell_basis.eval_gradients(msh, cl, qp.point(), 0, cell_degree);
            assert(num_cell_dofs == dphi.rows());
            assert(dphi.cols() == DIM);

            for(size_t i = 0; i < num_grad_dofs; i++){
               for(size_t j = i; j < num_grad_dofs; j++){
                  MG(i,j) += qp.weight() * mm_prod(gphi[i], gphi[j]);
               }

               for(size_t j = 0; j < num_cell_dofs; j++){
                  auto dphi_j = gphi[i];
                  for(size_t d = 0; d < dphi.cols(); d++) { dphi_j(d) = dphi(j,d);}
                  BG(i,j) += qp.weight() * mm_prod(gphi[i], dphi_j );
               }
            }
         }


         // lower part MG
         for(size_t i = 1; i < num_grad_dofs; i++)
            for(size_t j = 0; j < i; j++)
               MG(i,j) = MG(j,i);

            for (size_t face_i = 0; face_i < num_faces; face_i++)
            {
               auto current_face_range = dsr.face_range(face_i);
               auto fc = fcs[face_i];
               auto n = normal(msh, cl, fc);
               auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);

               for (auto& qp : face_quadpoints)
               {
                  matrix_type c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
                  auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());

                  // tau.n
                  matrix_type gphi_n = matrix_type::Zero(gphi.size(), 1);
                  for(size_t i = 0; i < gphi.size(); i++)
                     gphi_n(i,0) = mm_prod(gphi[i], n);


                  matrix_type T = qp.weight() * gphi_n * c_phi.transpose();

                  assert(T.rows() == num_grad_dofs);
                  assert(T.cols() == num_cell_dofs);

                  BG.block(0, 0, num_grad_dofs, num_cell_dofs) -= T;

                  matrix_type f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point(), 0, face_degree);
                  matrix_type F = qp.weight() * gphi_n * f_phi.transpose();

                  assert(F.rows() == num_grad_dofs);
                  assert(F.cols() == current_face_range.size());

                  BG.block(0, current_face_range.min(),
                           num_grad_dofs, current_face_range.size()) += F;
               }
            }

            m_oper  = MG.llt().solve(BG);    // GT
            if(compute_data)
               m_data  = BG.transpose() * m_oper;  // A

            assert(m_oper.rows() == num_grad_dofs);
            assert(m_oper.cols() == dsr.total_size());
      }
   };


   template<typename BQData>
   class gradient_reconstruction_elas_full_bq
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;

      typedef dynamic_matrix<scalar_type>         matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;

      const BQData&                               m_bqd;

      matrix_type     m_oper;
      matrix_type     m_data;

   public:

      matrix_type oper() const {return m_oper;}
      matrix_type data() const {return m_data;}

      gradient_reconstruction_elas_full_bq(const BQData& bqd) : m_bqd(bqd)
      {}

      void compute(const mesh_type& msh, const cell_type& cl, const bool compute_data = true)
      {
         const size_t DIM = msh.dimension;
         const size_t DIM2 = DIM * DIM;
         const size_t cell_degree = m_bqd.cell_degree();
         const size_t face_degree = m_bqd.face_degree();
         const size_t grad_degree = m_bqd.grad_degree();
         const size_t num_cell_dofs = (m_bqd.cell_basis.range(0, cell_degree)).size();
         const size_t num_face_dofs = (m_bqd.face_basis.range(0, face_degree)).size();
         const size_t num_grad_dofs = DIM * (m_bqd.grad_basis.range(0, grad_degree)).size();

         auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();

         timecounter tc;
         double t_base(0.0); double t_cons(0.0); double t_inv(0.0);

         assert(num_grad_dofs == DIM * DIM * binomial( grad_degree +DIM, grad_degree));

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         assert(dsr.total_size() == (num_cell_dofs + num_faces * num_face_dofs));

         matrix_type BG = matrix_type::Zero(num_grad_dofs, dsr.total_size());

         matrix_type MG = matrix_type::Zero(num_grad_dofs, num_grad_dofs);


         auto grad_quadpoints = m_bqd.grad_quadrature.integrate(msh, cl);
         for (auto& qp : grad_quadpoints)
         {
            tc.tic();
            auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());
            assert(num_grad_dofs == gphi.size());

            auto dphi = m_bqd.cell_basis.eval_gradients(msh, cl, qp.point());
            assert(num_cell_dofs == dphi.size());
            tc.toc();
            t_base += tc.to_double();

            tc.tic();

            size_t col = 0;
            for(size_t j = 0; j < num_grad_dofs; j += DIM2 ){
               for(size_t k = 0; k < DIM; k++ ){//depend de l'ordre des bases
                  for(size_t l = 0; l < DIM; l++ ){//depend de l'ordre des bases
                     for(size_t i = col; i < num_grad_dofs; i += DIM2){
                        MG(i,col) += qp.weight() * gphi[i](l,k) * gphi[col](l,k);
                     }
                     col++;
                  }
               }
            }

// Works but not optimized
//             for(size_t i = 0; i < num_grad_dofs; i++){
//                for(size_t j = i; j < num_grad_dofs; j++){
//                  MG(i,j) += qp.weight() * mm_prod(gphi[i], gphi[j]);


            size_t row = 0;
            for(size_t i = 0; i < num_grad_dofs; i += DIM2 ){
               for(size_t k = 0; k < DIM; k++ ){//depend de l'ordre des bases
                  for(size_t l = 0; l < DIM; l++ ){//depend de l'ordre des bases
                     for(size_t j = l; j < num_cell_dofs; j += DIM){
                        BG(row,j) += qp.weight() * gphi[row](l,k) * dphi[j](l,k);
                     }
                     row++;
                  }
               }
            }

            // Works but not optimized
//             for(size_t i = 0; i < num_grad_dofs; i++){
//                for(size_t j = ; j < num_cell_dofs; j++){
//                   BG(i,j) += qp.weight() * mm_prod(gphi[i], dphi[j]);
//                }
//             }

            tc.toc();
            t_cons += tc.to_double();
         }

         tc.tic();
         // lower part MG
         for(size_t i = 1; i <  num_grad_dofs; i++)
            for(size_t j = i; j < num_grad_dofs; j++)
               MG(i,j) = MG(j,i);

            tc.toc();
         t_cons += tc.to_double();

            for (size_t face_i = 0; face_i < num_faces; face_i++)
            {
               auto current_face_range = dsr.face_range(face_i);
               auto fc = fcs[face_i];
               auto n = normal(msh, cl, fc);
               auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);

               for (auto& qp : face_quadpoints)
               {
                  tc.tic();
                  auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point());
                  auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());
                  tc.toc();
                  t_base += tc.to_double();

                  // tau.n
                  tc.tic();
                  decltype(c_phi) gphi_n;

                  //vector zero
                  auto vzero(c_phi[0]);
                  vzero.setZero();

                  gphi_n.reserve(gphi.size());

//                   for(size_t i= 0; i < gphi.size(); i++){
//                      gphi_n.push_back(mm_prod(gphi[i] , n));
//                   }

                  for(size_t i = 0; i < num_grad_dofs; i += DIM2 ){
                     for(size_t k = 0; k < DIM; k++ ){//depend de l'ordre des bases
                        scalar_type val = gphi[i](0,0) * n(k);
                        for(size_t l = 0; l < DIM; l++ ){//depend de l'ordre des bases
                           auto gn(vzero);
                           gn(l) = val;
                           gphi_n.push_back(gn);
                        }
                     }
                  }

                  assert(gphi_n.size() == num_grad_dofs);

                  matrix_type T = matrix_type::Zero(num_grad_dofs, num_cell_dofs);

                  for(size_t j = 0; j < num_cell_dofs; j += DIM){
                     size_t col = j;
                     for(size_t k = 0; k < DIM; k++){
                        for(size_t i = 0; i < num_grad_dofs; i++){
                           T(i,col) = qp.weight() * gphi_n[i](k) * c_phi[col](k);
                        }
                        col++;
                     }
                  }

 // Works but not optimized
//                   for(size_t j = 0; j < num_cell_dofs; j++){
//                      for(size_t i = 0; i < num_grad_dofs; i++){
//                         T(i,j) = qp.weight() * mm_prod(gphi_n[i], c_phi[j]);
//                      }
//                   }

                  assert(T.rows() == num_grad_dofs);
                  assert(T.cols() == num_cell_dofs);

                  BG.block(0, 0, num_grad_dofs, num_cell_dofs) -= T;
                  tc.toc();
                  t_cons += tc.to_double();

                  tc.tic();
                  auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());
                  tc.toc();
                  t_base += tc.to_double();

                  tc.tic();
                  matrix_type F = matrix_type::Zero(BG.rows(), current_face_range.size());

                  for(size_t j = 0; j < current_face_range.size(); j += DIM){
                     size_t col = j;
                     for(size_t k = 0; k < DIM; k++){
                        for(size_t i = 0; i < num_grad_dofs; i++){
                           F(i,col) = qp.weight() * gphi_n[i](k) * f_phi[col](k);
                        }
                        col++;
                     }
                  }

                  // Works but not optimized
//                   for(size_t j = 0; j < current_face_range.size(); j++){
//                      for(size_t i = 0; i < num_grad_dofs; i++){
//                         F(i,j) = qp.weight() * mm_prod(gphi_n[i], f_phi[j]);
//                      }
//                   }

                  assert(F.rows() == num_grad_dofs);
                  assert(F.cols() == current_face_range.size());

                  BG.block(0, current_face_range.min(),
                           num_grad_dofs, current_face_range.size()) += F;
                   tc.toc();
                   t_cons += tc.to_double();
               }
            }

            tc.tic();
            m_oper  = MG.llt().solve(BG);    // GT
            tc.toc();
            t_inv += tc.to_double();
            if(compute_data)
               m_data  = BG.transpose() * m_oper;  // A

//             std::cout << "t_base = " << t_base << " sec" << std::endl;
//             std::cout << "t_cons = " << t_cons << " sec" << std::endl;
//             std::cout << "t_inv = " << t_inv << " sec" << std::endl;

            assert(m_oper.rows() == num_grad_dofs);
            assert(m_oper.cols() == dsr.total_size());
      }
   };


   template<typename BQData>
   class assembler_nl_bq
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;

      typedef dynamic_matrix<scalar_type>         matrix_type;


      typedef Eigen::Triplet<scalar_type>         triplet_type;

      std::vector<triplet_type>                   m_triplets;
      size_t                                      m_num_unknowns;

      const BQData&                               m_bqd;

   public:

      typedef Eigen::SparseMatrix<scalar_type>    sparse_matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;

      sparse_matrix_type      matrix;
      vector_type             rhs;

      assembler_nl_bq()                 = delete;

      assembler_nl_bq(const mesh_type& msh, const BQData& bqd)
      : m_bqd(bqd)
      {
         const size_t face_degree = m_bqd.face_degree();
         const size_t num_face_dofs = (m_bqd.face_basis.range(0, face_degree)).size();

         m_num_unknowns = num_face_dofs * (msh.faces_size() + msh.boundary_faces_size());
         matrix = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
         rhs = vector_type::Zero(m_num_unknowns);
      }

      template<typename LocalContrib>
      void
      assemble(const mesh_type& msh, const cell_type& cl, const LocalContrib& lc)
      {
         auto fcs = faces(msh, cl);
         const size_t face_degree = m_bqd.face_degree();
         const size_t face_basis_size = (m_bqd.face_basis.range(0, face_degree)).size();

         std::vector<size_t> l2g(fcs.size() * face_basis_size);
         for (size_t face_i = 0; face_i < fcs.size(); face_i++)
         {
            auto fc = fcs[face_i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
               throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            auto face_offset = face_id * face_basis_size;

            auto pos = face_i * face_basis_size;

            for (size_t i = 0; i < face_basis_size; i++)
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
         const size_t face_degree = m_bqd.face_degree();
         const size_t face_basis_size = (m_bqd.face_basis.range(0, face_degree)).size();

         size_t face_i = 0;
         for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
         {
            auto bfc = *itor;

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first)
               throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            auto face_offset = face_id * face_basis_size;
            auto face_offset_lagrange = (msh.faces_size() + face_i) * face_basis_size;

            auto fqd = m_bqd.face_quadrature.integrate(msh, bfc);

            matrix_type MFF     = matrix_type::Zero(face_basis_size, face_basis_size);
            vector_type rhs_f   = vector_type::Zero(face_basis_size);

            for (auto& qp : fqd)
            {
               auto f_phi = m_bqd.face_basis.eval_functions(msh, bfc, qp.point());
               assert(f_phi.size() == face_basis_size);
               for(size_t i = 0; i < face_basis_size; i++)
                  for(size_t j = i; j < face_basis_size; j++)
                     MFF(i,j) += qp.weight() * mm_prod(f_phi[i], f_phi[j]);


                  for(size_t i = 0; i < face_basis_size; i++)
                     rhs_f(i) += qp.weight() * mm_prod( f_phi[i],  bc(qp.point()));
            }

            //lower part
            for(size_t i = 1; i < face_basis_size; i++)
               for(size_t j = 0; j < i; j++)
                  MFF(i,j) = MFF(j,i);


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


   template<typename T, int DIM>
   static_vector<T,DIM>
   compute_gradient_matrix_pt(const dynamic_vector<T>& gradrec_coeff,
                              const std::vector<static_vector<T, DIM> >& base)
   {
      static_vector<T, DIM> ret = static_vector<T, DIM>::Zero();
      assert(gradrec_coeff.size() == base.size());

      for(size_t i = 0; i < base.size(); i++){
         ret += gradrec_coeff(i) * base[i];
      }

      return ret;
   }


   template<typename BQData>
   class projector_elas_bq
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;

      typedef dynamic_matrix<scalar_type>         matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;

      const BQData&                               m_bqd;

   public:

      projector_elas_bq(const BQData& bqd) : m_bqd(bqd)
      {}

      matrix_type cell_mm;
      matrix_type whole_mm;
      matrix_type grad_mm;

      template<typename Function>
      vector_type
      compute_cell(const mesh_type& msh, const cell_type& cl, const Function& f)
      {
         const size_t cell_degree = m_bqd.cell_degree();
         const size_t cell_basis_size = (m_bqd.cell_basis.range(0, cell_degree)).size();

         matrix_type mm = matrix_type::Zero(cell_basis_size, cell_basis_size);
         vector_type rhs = vector_type::Zero(cell_basis_size);

         auto cell_quadpoints = m_bqd.cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            auto phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point());

            for(size_t i = 0; i < cell_basis_size; i++)
               for(size_t j = i; j < cell_basis_size; j++)
                  mm(i,j)  += qp.weight() * mm_prod(phi[i], phi[j]);

               //lower part
               for (size_t i = 1; i < cell_basis_size; i++)
                  for (size_t j = 0; j < i; j++)
                     mm(i,j) = mm(j,i);

                  for(size_t i=0; i < cell_basis_size; i++){
                     rhs(i) += qp.weight() * mm_prod( f(qp.point()) , phi[i]);
                  }
         }

         cell_mm = mm;
         return mm.llt().solve(rhs);
      }


      template<typename Function>
      vector_type
      compute_whole(const mesh_type& msh, const cell_type& cl, const Function& f)
      {
         const size_t cell_degree = m_bqd.cell_degree();
         const size_t face_degree = m_bqd.face_degree();
         const size_t cell_basis_size = (m_bqd.cell_basis.range(0, cell_degree)).size();
         const size_t face_basis_size = m_bqd.face_basis.size();
         auto fcs = faces(msh, cl);

         compute_cell(msh, cl, f);

         vector_type ret = vector_type::Zero(cell_basis_size + fcs.size()*face_basis_size);
         whole_mm = matrix_type::Zero(cell_basis_size + fcs.size()*face_basis_size,
                                      cell_basis_size + fcs.size()*face_basis_size);

         ret.block(0, 0, cell_basis_size, 1) = compute_cell(msh, cl, f);
         whole_mm.block(0, 0, cell_basis_size, cell_basis_size) = cell_mm;

         size_t face_offset = cell_basis_size;
         for (auto& fc : fcs)
         {
            matrix_type mm = matrix_type::Zero(face_basis_size, face_basis_size);
            vector_type rhs = vector_type::Zero(face_basis_size);

            auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
               auto phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());
               assert(phi.size() == face_basis_size);

               for(size_t i=0; i < face_basis_size; i++)
                  for(size_t j=0; j < face_basis_size; j++)
                     mm(i,j)  += qp.weight() *mm_prod(phi[i], phi[j]);

                  for(size_t i=0; i < face_basis_size; i++)
                     rhs(i) += qp.weight() * mm_prod( f(qp.point()) , phi[i]);

            }

            ret.block(face_offset, 0, face_basis_size, 1) = mm.llt().solve(rhs);
            whole_mm.block(face_offset, face_offset, face_basis_size, face_basis_size) = mm;
            face_offset += face_basis_size;
         }


         return ret;
      }

      template<typename Function>
      vector_type
      compute_cell_grad(const mesh_type& msh, const cell_type& cl, const Function& f)
      {
         const size_t DIM = msh.dimension;
         const size_t DIM2 = DIM * DIM;
         const size_t grad_degree = m_bqd.grad_degree();
         const size_t grad_basis_size = DIM * (m_bqd.grad_basis.range(0, grad_degree)).size();

         matrix_type mm = matrix_type::Zero(grad_basis_size, grad_basis_size);
         vector_type rhs = vector_type::Zero(grad_basis_size);

         auto grad_quadpoints = m_bqd.grad_quadrature.integrate(msh, cl);
         for (auto& qp : grad_quadpoints)
         {
            auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());


            for(size_t j = 0; j < grad_basis_size; j += DIM2 ){
               size_t col = j;
               for(size_t k = 0; k < DIM; k++ ){//depend de l'ordre des bases
                  for(size_t l = 0; l < DIM; l++ ){//depend de l'ordre des bases
                     for(size_t i = col; i < grad_basis_size; i += DIM2){
                        mm(i,col) += qp.weight() * gphi[i](l,k) * gphi[col](l,k);
                     }
                     col++;
                  }
               }
            }

            //lower part
            for (size_t i = 1; i < grad_basis_size; i++)
               for (size_t j = i; j < grad_basis_size; j++)
                  mm(i,j) = mm(j,i);

               for(size_t i=0; i < grad_basis_size; i++){
                  rhs(i) += qp.weight() * mm_prod( f(qp.point()) , gphi[i]);
               }

         }

         grad_mm = mm;
         return mm.llt().solve(rhs);
      }

   };



} // namespace disk
