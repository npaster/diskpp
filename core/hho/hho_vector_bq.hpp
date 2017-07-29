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

//    template<typename Mesh,
//    template<typename, typename> class BasisFunction,
//    template<typename, typename> class BasisGradient,
//    template<typename, typename> class Quadrature>
//    class basis_quadrature_data_full /* this name really sucks */
//    {
//    public:
//       typedef Mesh                            mesh_type;
//       typedef typename mesh_type::cell        cell_type;
//       typedef typename mesh_type::face        face_type;
//
//       typedef BasisFunction<mesh_type, cell_type>         cell_basis_type;
//       typedef BasisFunction<mesh_type, face_type>         face_basis_type;
//       typedef BasisGradient<mesh_type, cell_type>         grad_basis_type;
//
//       typedef Quadrature<mesh_type, cell_type>    cell_quad_type;
//       typedef Quadrature<mesh_type, face_type>    face_quad_type;
//       typedef Quadrature<mesh_type, cell_type>    grad_quad_type;
//
//       cell_basis_type     cell_basis;
//       face_basis_type     face_basis;
//       cell_quad_type      cell_quadrature;
//       face_quad_type      face_quadrature;
//       grad_basis_type     grad_basis;
//       grad_quad_type      grad_quadrature;
//
//    private:
//       size_t  m_cell_degree, m_face_degree;
//
//       void init(void)
//       {
//          cell_basis          = cell_basis_type(m_cell_degree + 1);
//          face_basis          = face_basis_type(m_face_degree);
//          grad_basis          = grad_basis_type(m_cell_degree);
//          cell_quadrature     = cell_quad_type(2 * (m_cell_degree + 1));
//          face_quadrature     = face_quad_type(2 * m_face_degree);
//          grad_quadrature     = grad_quad_type(2 * m_cell_degree);
//       }
//
//    public:
//       basis_quadrature_data_full() : m_cell_degree(1), m_face_degree(1)
//       {
//          init();
//       }
//
//       basis_quadrature_data_full(size_t cell_degree, size_t face_degree)
//       {
//          if ( (cell_degree + 1 < face_degree) or (cell_degree > face_degree + 1) )
//             throw std::invalid_argument("Invalid cell degree");
//
//          m_cell_degree = cell_degree;
//          m_face_degree = face_degree;
//
//          init();
//       }
//
//       size_t cell_degree(void) const { return m_cell_degree; }
//       size_t face_degree(void) const { return m_face_degree; }
//    };



   template<typename BQData>
   class gradient_reconstruction_elas_bq
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;

      typedef typename BQData::cell_basis_type    cell_basis_type;
      typedef typename BQData::cell_quad_type     cell_quad_type;

      typedef dynamic_matrix<scalar_type>         matrix_type;

      const BQData&                               m_bqd;

      cell_basis_type     cell_basis;
      cell_quad_type      cell_quadrature;

   public:
      matrix_type     oper;
      matrix_type     data;

      gradient_reconstruction_elas_bq(const BQData& bqd) : m_bqd(bqd)
      {
         cell_basis          = cell_basis_type(bqd.cell_degree() + 1);
         cell_quadrature     = cell_quad_type(2 * (bqd.cell_degree() + 1));
      }

      void compute(const mesh_type& msh, const cell_type& cl)
      {
         const size_t cell_degree = m_bqd.cell_degree();
         const size_t face_degree = m_bqd.face_degree();
         const size_t cell_basis_size = (cell_basis.range(0, cell_degree + 1)).size();
         const size_t face_basis_size = m_bqd.face_basis.size();


         const size_t DIM= msh.dimension;
         const size_t dpk1 = DIM * binomial(cell_degree +1 + DIM, cell_degree +1);
         const size_t dpk0 = DIM * binomial( DIM, 0);
         const size_t dpk = DIM * binomial(cell_degree  + DIM, cell_degree );
         const size_t dpkf = DIM * binomial(face_degree  + DIM -1, face_degree );

         matrix_type stiff_mat = matrix_type::Zero(cell_basis_size, cell_basis_size);

         auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            auto dphi = cell_basis.eval_gradients(msh, cl, qp.point());
            assert(cell_basis_size == dphi.size());

            for(size_t i = 0; i < cell_basis_size; i++){
               for(size_t j = i; j < cell_basis_size; j++){
                  stiff_mat(i,j) += qp.weight() * mm_prod(dphi[i], dphi[j]);
               }
            }
         }

         // lower part
         for(size_t i = 1; i < cell_basis_size; i++)
            for(size_t j = 0; j < i; j++)
               stiff_mat(i,j) = stiff_mat(j,i);

         /* LHS: take basis functions derivatives from degree 1 to K+1 */
         auto MG_rowcol_range = cell_basis.range(1, cell_degree + 1);
         assert(MG_rowcol_range.from() == dpk0);
         assert(MG_rowcol_range.size() == (dpk1 - dpk0));
         matrix_type MG = take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

         /* RHS, volumetric part. */
         auto BG_row_range = cell_basis.range(1, cell_degree + 1);
         auto BG_col_range = cell_basis.range(0, cell_degree);

         assert(BG_row_range.from() == (dpk0));
         assert(BG_col_range.from() == 0) ;
         assert(BG_row_range.size() == (dpk1 - dpk0));
         assert(BG_col_range.size() == dpk) ;

         auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();

         const size_t num_cell_dofs = BG_col_range.size();

         assert(num_cell_dofs == dpk);

         dofspace_ranges dsr(num_cell_dofs, face_basis_size, num_faces);

         assert(dsr.total_size() == (num_cell_dofs + num_faces *face_basis_size));

         matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

         BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
         take(stiff_mat, BG_row_range, BG_col_range);

         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto current_face_range = dsr.face_range(face_i);
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);

            auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);

            auto cell_range = dsr.cell_range();

            assert(cell_range.min() == 0);
            assert(cell_range.max() == dpk);
            assert(cell_range.size() == dpk);

            for (auto& qp : face_quadpoints)
            {
               auto c_phi = cell_basis.eval_functions(msh, cl, qp.point()); // 0, m_degree);
               auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point()); // 1, m_degree+1);

               assert(c_phi.size() == dpk1);
               assert(c_dphi.size() == dpk1);

               decltype(c_phi) c_dphi_n;

               c_dphi_n.reserve(BG_row_range.to() - BG_row_range.from());

               assert(BG_row_range.from() == dpk0);
               assert(BG_row_range.to() == dpk1);

               for(size_t i=BG_row_range.from(); i< BG_row_range.to(); i++){
                  c_dphi_n.push_back(mm_prod(c_dphi[i] , n));
               }

               assert(c_dphi_n.size() == (dpk1 - dpk0));

               matrix_type  T= matrix_type::Zero(BG.rows(), BG_col_range.size());

               assert(c_dphi_n.size() == BG.rows());

               for(size_t i=0; i< BG.rows(); i++){
                  for(size_t j=0; j<BG_col_range.size(); j++){
                     T(i,j) = qp.weight() * mm_prod(c_dphi_n[i], c_phi[j]);
                  }
               }

               BG.block(0, 0, BG.rows(), BG_col_range.size()) -= T;

               auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());

               assert(f_phi.size() == dpkf);
               assert(current_face_range.size() == dpkf);

               matrix_type  F = matrix_type::Zero(BG.rows(), current_face_range.size());

               for(size_t i=0; i< BG.rows(); i++){
                  for(size_t j=0; j < current_face_range.size(); j++){
                     F(i,j) = qp.weight() * mm_prod(c_dphi_n[i], f_phi[j]);
                  }
               }

               BG.block(0, current_face_range.min(),
                        BG.rows(), current_face_range.size()) += F;
            }
         }

         assert(MG.rows() ==MG.cols());
         assert(MG.cols() == BG.rows());

         oper  = MG.ldlt().solve(BG);    // GT
         data  = BG.transpose() * oper;  // A
      }
   };



//    template<typename BQData>
//    class gradient_reconstruction_elas_full_bq
//    {
//       typedef typename BQData::mesh_type          mesh_type;
//       typedef typename mesh_type::scalar_type     scalar_type;
//       typedef typename mesh_type::cell            cell_type;
//
//       typedef dynamic_matrix<scalar_type>         matrix_type;
//
//       const BQData&                               m_bqd;
//
//    public:
//       matrix_type     oper;
//       matrix_type     data;
//
//       gradient_reconstruction_elas_full_bq(const BQData& bqd) : m_bqd(bqd)
//       {}
//
//       void compute(const mesh_type& msh, const cell_type& cl)
//       {
//          const size_t DIM= msh.dimension;
//          const size_t cell_degree = m_bqd.cell_degree();
//          const size_t face_degree = m_bqd.face_degree();
//          const size_t cell_basis_size = (m_bqd.cell_basis.range(0, cell_degree)).size();
//          const size_t grad_basis_size = DIM * (m_bqd.cell_basis.range(0, cell_degree)).size();
//          const size_t face_basis_size = m_bqd.face_basis.size();
//
//
//          const size_t dpk = DIM * binomial(cell_degree  + DIM, cell_degree );
//          const size_t gpk = DIM * DIM * binomial(cell_degree  + DIM, cell_degree );
//          const size_t dpkf = DIM * binomial(face_degree  + DIM -1, face_degree );
//
//
//          assert(cell_basis_size == dpk);
//          assert(grad_basis_size == gpk);
//          assert(face_basis_size == dpkf);
//
//          auto fcs = faces(msh, cl);
//          const size_t num_faces = fcs.size();
//
//          dofspace_ranges dsr(cell_basis_size, face_basis_size, num_faces);
//
//          assert(dsr.total_size() == (cell_basis_size + num_faces * face_basis_size));
//
//          matrix_type BG = matrix_type::Zero(grad_basis_size, dsr.total_size());
//
//          matrix_type MG = matrix_type::Zero(grad_basis_size, grad_basis_size);
//
//
//          // on sur integre legèrement
//          auto cell_quadpoints = m_bqd.grad_quadrature.integrate(msh, cl);
//          for (auto& qp : cell_quadpoints)
//          {
//             auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());
//             assert(gpk == gphi.size());
//
//             auto dphi = m_bqd.cell_basis.eval_gradients(msh, cl, qp.point());
//
//             for(size_t i = 0; i < grad_basis_size; i++){
//                for(size_t j = i; j < grad_basis_size; j++){
//                   MG(i,j) += qp.weight() * mm_prod(gphi[i], gphi[j]);
//                }
//
//
//                for(size_t j = 0; j < cell_basis_size; j++){
//                   BG(i,j) += qp.weight() * mm_prod(gphi[i], dphi[j]);
//                }
//             }
//          }
//
//          // lower part MG
//          for(size_t i = 1; i < grad_basis_size; i++)
//             for(size_t j = 0; j < i; j++)
//                MG(i,j) = MG(j,i);
//
//
//             for (size_t face_i = 0; face_i < num_faces; face_i++)
//             {
//                auto current_face_range = dsr.face_range(face_i);
//                auto fc = fcs[face_i];
//                auto n = normal(msh, cl, fc);
//
//                auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);
//
//                auto cell_range = dsr.cell_range();
//
//                assert(cell_range.min() == 0);
//                assert(cell_range.max() == dpk);
//                assert(cell_range.size() == dpk);
//
//                for (auto& qp : face_quadpoints)
//                {
//                   auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point()); // 0, m_degree);
//                   auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());
//
//                   // tau.n
//                   decltype(c_phi) gphi_n;
//
//                   gphi_n.reserve(gphi.size());
//
//                   for(size_t i= 0; i < grad_basis_size; i++){
//                      gphi_n.push_back(mm_prod(gphi[i] , n));
//                   }
//
//                   assert(gphi_n.size() == gpk);
//
//                   matrix_type T = matrix_type::Zero(BG.rows(), cell_basis_size);
//
//                   assert(gphi_n.size() == BG.rows());
//
//                   for(size_t i = 0; i < BG.rows(); i++){
//                      for(size_t j = 0; j < cell_basis_size; j++){
//                         T(i,j) = qp.weight() * mm_prod(gphi_n[i], c_phi[j]);
//                      }
//                   }
//
//                   BG.block(0, 0, BG.rows(), cell_basis_size) -= T;
//
//                   auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());
//
//                   assert(f_phi.size() == dpkf);
//                   assert(current_face_range.size() == dpkf);
//
//                   matrix_type  F= matrix_type::Zero(BG.rows(), current_face_range.size());
//
//                   for(size_t i=0; i< BG.rows(); i++){
//                      for(size_t j=0; j<current_face_range.size(); j++){
//                         F(i,j) = qp.weight() * mm_prod(gphi_n[i], f_phi[j]);
//                      }
//                   }
//
//                   BG.block(0, current_face_range.min(),
//                            BG.rows(), current_face_range.size()) += F;
//                }
//             }
//
//             assert(MG.rows() ==MG.cols());
//             assert(MG.cols() == BG.rows());
//
//             oper  = MG.ldlt().solve(BG);    // GT
//             data  = BG.transpose() * oper;  // A
//       }
//    };



   template<typename BQData>
   class elas_like_stabilization_bq
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;

      typedef dynamic_matrix<scalar_type>         matrix_type;

      typedef typename BQData::cell_basis_type    cell_basis_type;
      typedef typename BQData::cell_quad_type     cell_quad_type;

      const BQData&                               m_bqd;

      cell_basis_type     cell_basis;
      cell_quad_type      cell_quadrature;

   public:
      matrix_type     data;

      elas_like_stabilization_bq(const BQData& bqd) : m_bqd(bqd)
      {
         cell_basis          = cell_basis_type(bqd.cell_degree() + 1);
         cell_quadrature     = cell_quad_type(2 * (bqd.cell_degree() + 1));
      }

      void compute(const mesh_type& msh, const cell_type& cl, const matrix_type& gradrec_oper)
      {
         const size_t cell_degree = m_bqd.cell_degree();
         const size_t face_degree = m_bqd.face_degree();
         const size_t cell_basis_size = (cell_basis.range(0, cell_degree + 1)).size();
         const size_t face_basis_size = (m_bqd.face_basis.range(0, face_degree)).size();

         const size_t DIM= msh.dimension;
         const size_t dpk1 = DIM * binomial(cell_degree +1 + DIM, cell_degree +1);
         const size_t dpk0 = DIM * binomial( DIM, 0);
         const size_t dpk = DIM * binomial(cell_degree  + DIM, cell_degree );
         const size_t dpkf = DIM * binomial(face_degree  + DIM -1, face_degree );

         assert(cell_basis_size == dpk1);

         matrix_type mass_mat = matrix_type::Zero(cell_basis_size, cell_basis_size);

         auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
            assert(c_phi.size() == dpk1);

            for(size_t j = 0; j < cell_basis_size; j += DIM ){
               size_t col = j;
               for(size_t k = 0; k < DIM; k++ ){//depend de l'ordre des bases
                  for(size_t i = col; i < cell_basis_size; i += DIM){
                     mass_mat(i,col) += qp.weight() * c_phi[i](k) * c_phi[col](k);
                  }
                  col++;
               }
            }
         }

         //lower part
         for (size_t i = 0; i < cell_basis_size; i++)
            for (size_t j = i; j < cell_basis_size; j++)
               mass_mat(i,j) = mass_mat(j,i);

         auto zero_range         = cell_basis.range(0, cell_degree);
         auto one_range          = cell_basis.range(1, cell_degree+1);

         assert(zero_range.size() == dpk);
         assert(one_range.size() == (dpk1 - dpk0));

         // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

         //Step 1: compute \pi_T^k p_T^k v (third term).
         matrix_type M1 = take(mass_mat, zero_range, zero_range);
         matrix_type M2 = take(mass_mat, zero_range, one_range);
         matrix_type proj1 = -M1.llt().solve(M2*gradrec_oper);

         //Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
         matrix_type I_T = matrix_type::Identity(zero_range.size(), zero_range.size());
         proj1.block(0, 0, zero_range.size(), zero_range.size()) += I_T;

         auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();

         const size_t num_cell_dofs = (cell_basis.range(0, cell_degree)).size();

         dofspace_ranges dsr(num_cell_dofs, face_basis_size, num_faces);

         data = matrix_type::Zero(dsr.total_size(), dsr.total_size());

         // Step 3: project on faces (eqn. 21)
         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto current_face_range = dsr.face_range(face_i);
            auto h = diameter(msh, /*fcs[face_i]*/cl);
            auto fc = fcs[face_i];

            assert(face_basis_size == dpkf);

            matrix_type face_mass_matrix    = matrix_type::Zero(face_basis_size, face_basis_size);
            matrix_type face_trace_matrix   = matrix_type::Zero(face_basis_size, cell_basis_size);

            auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
               auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());
               auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
               //auto q_f_phi = qp.weight() * f_phi;

               assert(f_phi.size() == face_basis_size);
               assert(c_phi.size() == cell_basis_size);

               for(size_t j = 0; j < face_basis_size; j += DIM ){
                  size_t col = j;
                  for(size_t k = 0; k < DIM; k++ ){//depend de l'ordre des bases
                     for(size_t i = col; i < face_basis_size; i += DIM){
                        face_mass_matrix(i,col) += qp.weight() * f_phi[i](k) * f_phi[col](k);
                     }
                     col++;
                  }
               }

                  for(size_t i=0; i< face_basis_size; i++){
                     for(size_t j=0; j< cell_basis_size; j++){
                        face_trace_matrix(i,j) += qp.weight() * mm_prod(f_phi[i], c_phi[j]);
                     }
                  }
            }

            //lower part
            for (size_t i = 0; i < face_basis_size; i++)
               for (size_t j = i; j < face_basis_size; j++)
                  face_mass_matrix(i,j) = face_mass_matrix(j,i);

            Eigen::LLT<matrix_type> piKF;
            piKF.compute(face_mass_matrix);

            // Step 3a: \pi_F^k( v_F - p_T^k v )
            auto face_range = current_face_range.remove_offset();
            matrix_type MR1 = take(face_trace_matrix, face_range, one_range);
            assert(MR1.cols() == gradrec_oper.rows());
            matrix_type proj2 = piKF.solve(MR1*gradrec_oper);
            matrix_type I_F = matrix_type::Identity(face_basis_size, face_basis_size);
            auto block_offset = current_face_range.min();
            proj2.block(0, block_offset, face_basis_size, face_basis_size) -= I_F;

            // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
            matrix_type MR2 = take(face_trace_matrix, face_range, zero_range);
            assert(MR2.cols() == proj1.rows());
            matrix_type proj3 = piKF.solve(MR2*proj1);

            assert(proj2.rows() == proj3.rows());
            assert(proj2.cols() == proj3.cols());

            matrix_type BRF = proj2 + proj3;

            data += BRF.transpose() * face_mass_matrix * BRF / h;
         }
      }
   };



   template<typename BQData>
   class elas_like_stabilization_PIKF_bq
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;

      typedef dynamic_matrix<scalar_type>         matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;

      const BQData&                               m_bqd;

   public:
      matrix_type     data;

      elas_like_stabilization_PIKF_bq(const BQData& bqd) : m_bqd(bqd)
      {}


      void compute(const mesh_type& msh, const cell_type& cl)
      {
         const size_t cell_degree = m_bqd.cell_degree();
         const size_t face_degree = m_bqd.face_degree();
         const size_t cell_basis_size = (m_bqd.cell_basis.range(0, cell_degree)).size();
         const size_t face_basis_size = m_bqd.face_basis.size();

         const size_t DIM = msh.dimension;
         const size_t cpk = DIM * binomial(cell_degree  + DIM, cell_degree);
         const size_t fpk = DIM * binomial(face_degree  + DIM -1, face_degree);

         assert(cpk == cell_basis_size);
         assert(fpk == face_basis_size);

         auto zero_range = m_bqd.cell_basis.range(0, cell_degree);

         auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();

         dofspace_ranges dsr(cell_basis_size, face_basis_size, num_faces);

         data = matrix_type::Zero(dsr.total_size(), dsr.total_size());

         //Step 2: v_T
         matrix_type I_T = matrix_type::Identity(cell_basis_size, cell_basis_size);
         matrix_type proj1  = matrix_type::Zero(cell_basis_size, dsr.total_size());
         proj1.block(0, 0, cell_basis_size, cell_basis_size) += I_T;

         // Step 3: project on faces (eqn. 21)
         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto current_face_range = dsr.face_range(face_i);
            auto h = diameter(msh, /*fcs[face_i]*/cl);
            auto fc = fcs[face_i];

            matrix_type face_mass_matrix    = matrix_type::Zero(face_basis_size, face_basis_size);
            matrix_type face_trace_matrix   = matrix_type::Zero(face_basis_size, cell_basis_size);

            auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
               auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());
               auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point());

               assert(f_phi.size() == face_basis_size);

               for(size_t i = 0; i < face_basis_size; i++){
                  for(size_t j = i; j < face_basis_size; j++){
                     face_mass_matrix(i,j) += qp.weight() * mm_prod(f_phi[i], f_phi[j]);
                  }
               }

               //lower part
               for (size_t i = 1; i < face_basis_size; i++)
                  for (size_t j = 0; j < i; j++)
                     face_mass_matrix(i,j) = face_mass_matrix(j,i);

                  for(size_t i=0; i< face_basis_size; i++){
                     for(size_t j=0; j< cell_basis_size; j++){
                        face_trace_matrix(i,j) += qp.weight() * mm_prod(f_phi[i], c_phi[j]);
                     }
                  }
            }

            Eigen::LLT<matrix_type> piKF;
            piKF.compute(face_mass_matrix);

            // Step 3a: v_F
            auto face_range = current_face_range.remove_offset();

            matrix_type I_F = matrix_type::Identity(face_basis_size, face_basis_size);
            auto block_offset = current_face_range.min();
            matrix_type proj2  = matrix_type::Zero(face_basis_size, dsr.total_size());
            proj2.block(0, block_offset, face_basis_size, face_basis_size) -= I_F;

            // Step 3b: \pi_F^k( v_T )
            matrix_type MR2 = take(face_trace_matrix, face_range, zero_range);
            assert(MR2.cols() == proj1.rows());
            matrix_type proj3 = piKF.solve(MR2*proj1);

            assert(proj2.rows() == proj3.rows());
            assert(proj2.cols() == proj3.cols());

            // \pi_F^k( v_T ) - v_F
            matrix_type BRF = proj2 + proj3;

            data += BRF.transpose() * face_mass_matrix * BRF / h;
         }
      }
   };


   template<typename BQData>
   class elas_like_stabilization_L2_bq
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;

      typedef dynamic_matrix<scalar_type>         matrix_type;

      typedef typename BQData::cell_basis_type    cell_basis_type;
      typedef typename BQData::cell_quad_type     cell_quad_type;

      const BQData&                               m_bqd;

      cell_basis_type     cell_basis;
      cell_quad_type      cell_quadrature;

   public:
      matrix_type     data;

      elas_like_stabilization_L2_bq(const BQData& bqd) : m_bqd(bqd)
      {
         cell_basis          = cell_basis_type(bqd.cell_degree() + 1);
         cell_quadrature     = cell_quad_type(2 * (bqd.cell_degree() + 1));
      }

      void compute(const mesh_type& msh, const cell_type& cl)
      {
         const size_t cell_degree = m_bqd.cell_degree();
         const size_t face_degree = m_bqd.face_degree();
         const size_t cell_basis_size = (m_bqd.cell_basis.range(0, cell_degree)).size();
         const size_t face_basis_size = m_bqd.face_basis.size();

         const size_t DIM = msh.dimension;
         const size_t cpk = DIM * binomial(cell_degree  + DIM, cell_degree);
         const size_t fpk = DIM * binomial(face_degree  + DIM -1, face_degree);

         assert(cpk == cell_basis_size);
         assert(fpk == face_basis_size);

         auto zero_range = m_bqd.cell_basis.range(0, cell_degree);

         auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();

         dofspace_ranges dsr(cell_basis_size, face_basis_size, num_faces);

         data = matrix_type::Zero(dsr.total_size(), dsr.total_size());

         //Step 1: v_T
         matrix_type I_T = matrix_type::Identity(cell_basis_size, cell_basis_size);
         matrix_type proj1  = matrix_type::Zero(cell_basis_size, dsr.total_size());
         proj1.block(0, 0, cell_basis_size, cell_basis_size) += I_T;

         // Step 3: project on faces (eqn. 21)
         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto current_face_range = dsr.face_range(face_i);
            auto h = diameter(msh, /*fcs[face_i]*/cl);
            auto fc = fcs[face_i];

            matrix_type face_mass_matrix    = matrix_type::Zero(face_basis_size, face_basis_size);
            matrix_type face_trace_matrix   = matrix_type::Zero(face_basis_size, cell_basis_size);

            auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
               auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());
               auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point());

               assert(f_phi.size() == face_basis_size);

               for(size_t i = 0; i < face_basis_size; i++){
                  for(size_t j = i; j < face_basis_size; j++){
                     face_mass_matrix(i,j) += qp.weight() * mm_prod(f_phi[i], f_phi[j]);
                  }
               }

               //lower part
               for (size_t i = 1; i < face_basis_size; i++)
                  for (size_t j = 0; j < i; j++)
                     face_mass_matrix(i,j) = face_mass_matrix(j,i);

                  for(size_t i=0; i< face_basis_size; i++){
                     for(size_t j=0; j< cell_basis_size; j++){
                        face_trace_matrix(i,j) += qp.weight() * mm_prod(f_phi[i], c_phi[j]);
                     }
                  }
            }

            Eigen::LLT<matrix_type> piKF;
            piKF.compute(face_mass_matrix);

            // Step 3a: v_F
            auto face_range = current_face_range.remove_offset();

            matrix_type I_F = matrix_type::Identity(face_basis_size, face_basis_size);
            auto block_offset = current_face_range.min();
            matrix_type proj2  = matrix_type::Zero(face_basis_size, dsr.total_size());
            proj2.block(0, block_offset, face_basis_size, face_basis_size) -= I_F;

            // Step 3b: \pi_F^k( v_T )
            matrix_type MR2 = take(face_trace_matrix, face_range, zero_range);
            assert(MR2.cols() == proj1.rows());
            matrix_type proj3 = MR2*proj1;

            assert(proj2.rows() == proj3.rows());
            assert(proj2.cols() == proj3.cols());

            // \pi_F^k( v_T ) - v_F
            matrix_type BRF = proj2 + proj3;

            data += BRF.transpose() * face_mass_matrix * BRF;
         }
      }
   };


} // namespace disk
