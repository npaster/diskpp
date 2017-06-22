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


   template<typename Mesh, typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType>
   class sgradient_reconstruction_elas
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

      sgradient_reconstruction_elas()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }

      sgradient_reconstruction_elas(size_t degree)
      : m_degree(degree)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }


      void compute(const mesh_type& msh, const cell_type& cl)
      {
         const size_t DIM= msh.dimension;
         const size_t dpk1 = DIM * binomial(m_degree +1 + DIM, m_degree +1);
         const size_t dpk0 = DIM * binomial( DIM, 0);
         const size_t dpk = DIM * binomial(m_degree  + DIM, m_degree );
         const size_t dpkf = DIM * binomial(m_degree  + DIM -1, m_degree );
         
         //Maybe we can improve it
         auto SS = cell_basis.eval_ssgradients_const(msh, cl);
         
         const size_t offset_lag = cell_basis.size();
         const size_t nb_lag = SS.size();
         
         matrix_type stiff_mat = matrix_type::Zero(cell_basis.size() + nb_lag, cell_basis.size() + nb_lag);
         

         auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            auto sdphi = cell_basis.eval_sgradients(msh, cl, qp.point());
            auto dphi = cell_basis.eval_gradients(msh, cl, qp.point());
            assert(cell_basis.size() == sdphi.size());
            assert(dphi.size() == sdphi.size());
            assert(dpk1 == sdphi.size());
            for(size_t i = 0; i < sdphi.size(); i++){
               for(size_t j = i; j < sdphi.size(); j++){
                  stiff_mat(i,j) += qp.weight() * mm_prod(sdphi.at(i), sdphi.at(j));
               }
               
               // lagrangian part
               for(size_t j = 0; j < nb_lag; j++){
                  stiff_mat(i, offset_lag + j) +=  qp.weight() * mm_prod(dphi.at(i), SS.at(j));
               }
            }
         }

         
         
         // lower part
         for(size_t i = 1; i < stiff_mat.rows(); i++)
            for(size_t j = 0; j < i; j++)
               stiff_mat(i,j) = stiff_mat(j,i);

         /* LHS: take basis functions derivatives from degree 1 to K+1 */
         auto MG_rowcol_range = cell_basis.range(1, m_degree+1);
         assert(MG_rowcol_range.from() == dpk0);
         assert(MG_rowcol_range.size() == (dpk1 - dpk0));
         
         dof_range MG_rowcol_range_lag(MG_rowcol_range.from(), MG_rowcol_range.to() + nb_lag);
         matrix_type MG = take(stiff_mat, MG_rowcol_range_lag, MG_rowcol_range_lag);
         
         assert(MG.rows() == (dpk1 - dpk0) + nb_lag);
         assert(MG.rows() == MG.cols());

         /* RHS, volumetric part. */
         auto BG_row_range = cell_basis.range(1, m_degree+1);
         auto BG_col_range = cell_basis.range(0, m_degree);

         assert(BG_row_range.from() == (dpk0));
         assert(BG_col_range.from() == 0) ;
         assert(BG_row_range.size() == (dpk1 - dpk0));
         assert(BG_col_range.size() == dpk) ;

         auto fcs = faces(msh, cl);
         auto num_faces = fcs.size();

         auto num_cell_dofs = cell_basis.range(0, m_degree).size();
         auto num_face_dofs = face_basis.size();

         assert(num_cell_dofs == dpk);
         assert(num_face_dofs == dpkf);

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         assert(dsr.total_size() == (num_cell_dofs + num_faces * num_face_dofs));

         matrix_type BG = matrix_type::Zero(BG_row_range.size() + nb_lag, dsr.total_size());

         BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
         take(stiff_mat, BG_row_range, BG_col_range);

         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto current_face_range = dsr.face_range(face_i);
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);

            auto face_quadpoints = face_quadrature.integrate(msh, fc);

            auto cell_range = dsr.cell_range();

            assert(cell_range.min() == 0);
            assert(cell_range.max() == dpk);
            assert(cell_range.size() == dpk);

            for (auto& qp : face_quadpoints)
            {
               auto c_phi = cell_basis.eval_functions(msh, cl, qp.point()); // 0, m_degree);
               auto c_sdphi = cell_basis.eval_sgradients(msh, cl, qp.point()); // 1, m_degree+1);

               assert(c_phi.size() == dpk1);
               assert(c_sdphi.size() == dpk1);

               decltype(c_phi) c_sdphi_n;
               
               c_sdphi_n.reserve(BG_row_range.to() - BG_row_range.from());

               assert(BG_row_range.from() == dpk0);
               assert(BG_row_range.to() == dpk1);

               for(size_t i=BG_row_range.from(); i< BG_row_range.to(); i++){
                  c_sdphi_n.push_back(mm_prod(c_sdphi.at(i) , n));
               }

               assert(c_sdphi_n.size() == (dpk1 - dpk0));

               matrix_type  T= matrix_type::Zero(BG.rows(), BG_col_range.size());

               assert((c_sdphi_n.size() + nb_lag) == BG.rows());

               for(size_t i=0; i< c_sdphi_n.size(); i++){
                  for(size_t j=0; j < BG_col_range.size(); j++){
                     T(i,j) = qp.weight() * mm_prod(c_sdphi_n.at(i), c_phi.at(j));
                  }
               }

               BG.block(0, 0, BG.rows(), BG_col_range.size()) -= T;

               auto f_phi = face_basis.eval_functions(msh, fc, qp.point());

               assert(f_phi.size() == dpkf);
               assert(current_face_range.size() == dpkf);

               matrix_type  F= matrix_type::Zero(BG.rows(), current_face_range.size());

               for(size_t i=0; i< c_sdphi_n.size(); i++){
                  for(size_t j=0; j < current_face_range.size(); j++){
                     F(i,j) = qp.weight() * mm_prod(c_sdphi_n.at(i), f_phi.at(j));
                  }
               }

               BG.block(0, current_face_range.min(),
                        BG.rows(), current_face_range.size()) += F;
            }
         }

         assert(MG.rows() == MG.cols());
         assert(MG.cols() == BG.rows());
         
         oper  = (MG.ldlt().solve(BG)).topLeftCorner(MG_rowcol_range.size(), dsr.total_size());    // GT
         data  = (BG.topLeftCorner(MG_rowcol_range.size(), dsr.total_size())).transpose() * oper;  // A
      }
   };


   template<typename Mesh, typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType>
   class gradient_reconstruction_elas
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



      typedef material_tensor<scalar_type, mesh_type::dimension, mesh_type::dimension>
      material_tensor_type;

      cell_basis_type                             cell_basis;
      cell_quadrature_type                        cell_quadrature;

      face_basis_type                             face_basis;
      face_quadrature_type                        face_quadrature;

      size_t                                      m_degree;

   public:
      matrix_type     oper;
      matrix_type     data;

      gradient_reconstruction_elas()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }

      gradient_reconstruction_elas(size_t degree)
      : m_degree(degree)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }

      void compute(const mesh_type& msh, const cell_type& cl)
      {
         material_tensor_type id_tens;
         id_tens = material_tensor_type::Identity();
         compute(msh, cl, id_tens);
      }

      void compute(const mesh_type& msh, const cell_type& cl,
                   const material_tensor_type& mtens)
      {
         const size_t DIM= msh.dimension;
         const size_t dpk1 = DIM * binomial(m_degree +1 + DIM, m_degree +1);
         const size_t dpk0 = DIM * binomial( DIM, 0);
         const size_t dpk = DIM * binomial(m_degree  + DIM, m_degree );
         const size_t dpkf = DIM * binomial(m_degree  + DIM -1, m_degree );

         matrix_type stiff_mat = matrix_type::Zero(cell_basis.size(), cell_basis.size());

         auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            auto dphi = cell_basis.eval_gradients(msh, cl, qp.point());
            assert(cell_basis.size() == dphi.size());
            assert(dpk1 == dphi.size());
            for(size_t i = 0; i < dphi.size(); i++){
               for(size_t j = i; j < dphi.size(); j++){
                  stiff_mat(i,j) += qp.weight() * mm_prod(dphi[i], dphi[j]);
               }
            }
         }

         // lower part
         for(size_t i = 1; i < cell_basis.size(); i++)
            for(size_t j = 0; j < i; j++)
               stiff_mat(i,j) = stiff_mat(j,i);

            /* LHS: take basis functions derivatives from degree 1 to K+1 */
            auto MG_rowcol_range = cell_basis.range(1, m_degree+1);
         assert(MG_rowcol_range.from() == dpk0);
         assert(MG_rowcol_range.size() == (dpk1 - dpk0));
         matrix_type MG = take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

         /* RHS, volumetric part. */
         auto BG_row_range = cell_basis.range(1, m_degree+1);
         auto BG_col_range = cell_basis.range(0, m_degree);

         assert(BG_row_range.from() == (dpk0));
         assert(BG_col_range.from() == 0) ;
         assert(BG_row_range.size() == (dpk1 - dpk0));
         assert(BG_col_range.size() == dpk) ;

         auto fcs = faces(msh, cl);
         auto num_faces = fcs.size();

         auto num_cell_dofs = cell_basis.range(0, m_degree).size();
         auto num_face_dofs = face_basis.size();

         assert(num_cell_dofs == dpk);
         assert(num_face_dofs == dpkf);

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         assert(dsr.total_size() == (num_cell_dofs + num_faces * num_face_dofs));

         matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

         BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
         take(stiff_mat, BG_row_range, BG_col_range);

         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto current_face_range = dsr.face_range(face_i);
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);

            auto face_quadpoints = face_quadrature.integrate(msh, fc);

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

               auto f_phi = face_basis.eval_functions(msh, fc, qp.point());

               assert(f_phi.size() == dpkf);
               assert(current_face_range.size() == dpkf);

               matrix_type  F= matrix_type::Zero(BG.rows(), current_face_range.size());

               for(size_t i=0; i< BG.rows(); i++){
                  for(size_t j=0; j<current_face_range.size(); j++){
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

   
   template<typename Mesh, typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType,
   typename TestBasisType, typename TestQuadType>
   class gradient_reconstruction_elas_full
   {
      typedef Mesh                                mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;
      typedef typename mesh_type::face            face_type;
      
      typedef CellBasisType                       cell_basis_type;
      typedef CellQuadType                        cell_quadrature_type;
      typedef FaceBasisType                       face_basis_type;
      typedef FaceQuadType                        face_quadrature_type;
      typedef TestBasisType                       test_basis_type;
      typedef TestQuadType                        test_quadrature_type;
      
      typedef dynamic_matrix<scalar_type>         matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;
      
      
      cell_basis_type                             cell_basis;
      cell_quadrature_type                        cell_quadrature;
      
      face_basis_type                             face_basis;
      face_quadrature_type                        face_quadrature;
      
      test_basis_type                             test_basis;
      test_quadrature_type                        test_quadrature;
      
      size_t                                      m_degree;
      
   public:
      matrix_type     oper;
      matrix_type     data;
      
      gradient_reconstruction_elas_full()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree);
         cell_quadrature     = cell_quadrature_type(2*m_degree);
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
         test_basis          = test_basis_type(m_degree);
         test_quadrature     = test_quadrature_type(2*m_degree);
      }
      
      gradient_reconstruction_elas_full(size_t degree)
      : m_degree(degree)
      {
         cell_basis          = cell_basis_type(m_degree);
         cell_quadrature     = cell_quadrature_type(2*m_degree+1);
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
         test_basis          = test_basis_type(m_degree);
         test_quadrature     = test_quadrature_type(2*m_degree);
      }
      
      
      void compute(const mesh_type& msh, const cell_type& cl)
      {
         const size_t DIM= msh.dimension;
         const size_t dpk = DIM * binomial(m_degree  + DIM, m_degree );
         const size_t tpk = DIM * DIM * binomial(m_degree  + DIM, m_degree );
         const size_t dpkf = DIM * binomial(m_degree  + DIM -1, m_degree );
         
         //auto BG_row_range = test_basis.range(0, m_degree);
         auto BG_col_range = cell_basis.range(0, m_degree);
         
         assert(BG_col_range.size() == dpk);
         
         auto fcs = faces(msh, cl);
         auto num_faces = fcs.size();
         
         auto num_cell_dofs = cell_basis.range(0, m_degree).size();
         auto num_face_dofs = face_basis.size();
         
         assert(num_cell_dofs == dpk);
         assert(num_face_dofs == dpkf);
         
         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);
         
         assert(dsr.total_size() == (num_cell_dofs + num_faces * num_face_dofs));
         
         matrix_type BG = matrix_type::Zero(tpk, dsr.total_size());
         
         matrix_type MG = matrix_type::Zero(tpk, tpk);
         
         
         // on sur integre legèrement
         auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            auto tphi = test_basis.eval_functions(msh, cl, qp.point());
            assert(tpk == tphi.size());
            
            auto dphi = cell_basis.eval_gradients(msh, cl, qp.point());
            assert(cell_basis.size() == dphi.size());
            assert(dpk == dphi.size());
            for(size_t i = 0; i < tphi.size(); i++){
               for(size_t j = i; j < tphi.size(); j++){
                  MG(i,j) += qp.weight() * mm_prod(tphi[i], tphi[j]);
               }
               
               
               for(size_t j = 0; j < dphi.size(); j++){
                  BG(i,j) += qp.weight() * mm_prod(tphi[i], dphi[j]);
               }
            }
         }

         // lower part MG
         for(size_t i = 1; i < tpk; i++)
            for(size_t j = 0; j < i; j++)
               MG(i,j) = MG(j,i);
 
         
         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto current_face_range = dsr.face_range(face_i);
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);
            
            auto face_quadpoints = face_quadrature.integrate(msh, fc);
            
            auto cell_range = dsr.cell_range();
            
            assert(cell_range.min() == 0);
            assert(cell_range.max() == dpk);
            assert(cell_range.size() == dpk);
            
            for (auto& qp : face_quadpoints)
            {
               auto c_phi = cell_basis.eval_functions(msh, cl, qp.point()); // 0, m_degree);
               auto tphi = test_basis.eval_functions(msh, cl, qp.point());
               
               assert(c_phi.size() == dpk);
               
               // tau.n
               decltype(c_phi) tphi_n;
               
               tphi_n.reserve(tphi.size());
               
               for(size_t i= 0; i < tphi.size(); i++){
                  tphi_n.push_back(mm_prod(tphi[i] , n));
               }
               
               assert(tphi_n.size() == tpk);
               
               matrix_type  T= matrix_type::Zero(BG.rows(), BG_col_range.size());
               
               assert(tphi_n.size() == BG.rows());
               
               for(size_t i=0; i< BG.rows(); i++){
                  for(size_t j=0; j<BG_col_range.size(); j++){
                     T(i,j) = qp.weight() * mm_prod(tphi_n[i], c_phi[j]);
                  }
               }
               
               BG.block(0, 0, BG.rows(), BG_col_range.size()) -= T;
               
               auto f_phi = face_basis.eval_functions(msh, fc, qp.point());
               
               assert(f_phi.size() == dpkf);
               assert(current_face_range.size() == dpkf);
               
               matrix_type  F= matrix_type::Zero(BG.rows(), current_face_range.size());
               
               for(size_t i=0; i< BG.rows(); i++){
                  for(size_t j=0; j<current_face_range.size(); j++){
                     F(i,j) = qp.weight() * mm_prod(tphi_n[i], f_phi[j]);
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


   template<typename Mesh, typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType,
   typename DivCellBasisType, typename DivCellQuadType>
   class divergence_reconstruction_elas
   {
      typedef Mesh                                mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;
      typedef typename mesh_type::face            face_type;

      typedef CellBasisType                       cell_basis_type;
      typedef CellQuadType                        cell_quadrature_type;
      typedef FaceBasisType                       face_basis_type;
      typedef FaceQuadType                        face_quadrature_type;
      typedef DivCellBasisType                    div_cell_basis_type;
      typedef DivCellQuadType                     div_cell_quadrature_type;

      typedef dynamic_matrix<scalar_type>         matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;

      cell_basis_type                             cell_basis;
      cell_quadrature_type                        cell_quadrature;

      face_basis_type                             face_basis;
      face_quadrature_type                        face_quadrature;

      div_cell_basis_type                         div_cell_basis;
      div_cell_quadrature_type                    div_cell_quadrature;

      size_t                                      m_degree;

   public:
      matrix_type     oper;
      matrix_type     data;

      divergence_reconstruction_elas()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
         div_cell_basis      = div_cell_basis_type(m_degree);
         div_cell_quadrature = div_cell_quadrature_type(2*m_degree);
      }

      divergence_reconstruction_elas(size_t degree)
      : m_degree(degree)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
         div_cell_basis      = div_cell_basis_type(m_degree);
         div_cell_quadrature = div_cell_quadrature_type(2*m_degree);
      }

      void compute(const mesh_type& msh, const cell_type& cl)
      {
         const size_t DIM= msh.dimension;
         const size_t dpk = DIM * binomial(m_degree  + DIM, m_degree );
         const size_t dpkf = DIM * binomial(m_degree  + DIM -1, m_degree );
         const size_t dpkd = binomial(m_degree  + DIM, m_degree );

         auto dcbs = div_cell_basis.size();
         matrix_type MD = matrix_type::Zero(dcbs, dcbs);

         assert(dcbs == dpkd);

         auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            auto phi_d = div_cell_basis.eval_functions(msh, cl, qp.point());
            assert(phi_d.size() == dpkd);
            for (size_t i = 0; i < dcbs; i++)
               for (size_t j = i; j < dcbs; j++)
                  MD(i,j) += qp.weight() * mm_prod(phi_d(i), phi_d(j));
         }

         //lower part
         for (size_t i = 1; i < dcbs; i++)
            for (size_t j = 0; j < i; j++)
               MD(i,j) = MD(j,i);

         /* RHS, volumetric part. */
         auto fcs = faces(msh, cl);
         auto num_cell_dofs = cell_basis.range(0, m_degree).size();
         auto num_face_dofs = face_basis.size();
         auto num_faces = fcs.size();

         assert(num_cell_dofs == dpk);

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         assert(dsr.total_size() == (num_cell_dofs + num_faces * num_face_dofs));

         matrix_type BD = matrix_type::Zero(dcbs, dsr.total_size());
         for (auto& qp : cell_quadpoints)
         {
            auto phi = cell_basis.eval_functions(msh, cl, qp.point());
            auto dphi_d = div_cell_basis.eval_gradients(msh, cl, qp.point());

            assert(dcbs == dphi_d.rows());
            assert(DIM == dphi_d.cols());


            for (size_t i = 0; i < dcbs; i++){
               for (size_t j = 0; j < num_cell_dofs; j++){
                  auto dphi_di(phi.at(j));
                  dphi_di.setZero();
                  assert(dphi_d.cols() == DIM);
                  for(size_t id=0; id< DIM; id++){
                     dphi_di(id)=dphi_d(i,id);
                  }
                  BD(i,j) -= qp.weight() * mm_prod(dphi_di, phi.at(j));
               }
            }
         }

         size_t face_offset = num_cell_dofs;

         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);
            auto face_quadpoints = face_quadrature.integrate(msh, fc);

            for (auto& qp : face_quadpoints)
            {
               auto phi_d = div_cell_basis.eval_functions(msh, cl, qp.point());
               auto phi = face_basis.eval_functions(msh, fc, qp.point());

               assert(phi_d.size() == dcbs);
               assert(phi.size() == dpkf);

               for (size_t i = 0; i < dcbs; i++) {
                  for (size_t j = 0; j < face_basis.size(); j++) {
                     scalar_type p1 = mm_prod(phi.at(j), n);
                     scalar_type p2 = p1 * phi_d(i);
                     BD(i,face_offset+j) += qp.weight() * p2;
                  }
               }
            }

            assert(face_basis.size() == dpkf);
            face_offset += face_basis.size();
         }
         oper = MD.ldlt().solve(BD);
         //oper = MD.partialPivLu().solve(BD);
         data = BD.transpose() * oper;

         assert(oper.rows() == dcbs);
         assert(data.rows() == dsr.total_size());
         assert(data.cols() == dsr.total_size());
      }
   };




   template<typename Mesh, typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType>
   class elas_like_stabilization
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

      elas_like_stabilization()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }

      elas_like_stabilization(size_t degree)
      : m_degree(degree)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }

      void compute(const mesh_type& msh, const cell_type& cl, const matrix_type& gradrec_oper)
      {
         const size_t DIM= msh.dimension;
         const size_t dpk1 = DIM * binomial(m_degree +1 + DIM, m_degree +1);
         const size_t dpk0 = DIM * binomial( DIM, 0);
         const size_t dpk = DIM * binomial(m_degree  + DIM, m_degree );
         const size_t dpkf = DIM * binomial(m_degree  + DIM -1, m_degree );

         assert(cell_basis.size() == dpk1);

         matrix_type mass_mat = matrix_type::Zero(cell_basis.size(), cell_basis.size());

         auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
            assert(c_phi.size() == dpk1);
            for(size_t i = 0; i < cell_basis.size(); i++){
               for(size_t j = i; j< cell_basis.size(); j++){
                  mass_mat(i,j) += qp.weight() * mm_prod(c_phi[i], c_phi[j]);
               }
            }
         }
         
         //lower part
         for (size_t i = 1; i < cell_basis.size(); i++)
            for (size_t j = 0; j < i; j++)
               mass_mat(i,j) = mass_mat(j,i);

         auto zero_range         = cell_basis.range(0, m_degree);
         auto one_range          = cell_basis.range(1, m_degree+1);

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
         auto num_faces = fcs.size();

         auto num_cell_dofs = cell_basis.range(0, m_degree).size();
         auto num_face_dofs = face_basis.size();

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         data = matrix_type::Zero(dsr.total_size(), dsr.total_size());

         // Step 3: project on faces (eqn. 21)
         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto current_face_range = dsr.face_range(face_i);
            auto fbs = current_face_range.size();
            auto h = diameter(msh, /*fcs[face_i]*/cl);
            auto fc = fcs[face_i];

            assert(fbs == dpkf);

            matrix_type face_mass_matrix    = matrix_type::Zero(fbs, fbs);
            matrix_type face_trace_matrix   = matrix_type::Zero(fbs, cell_basis.size());

            auto face_quadpoints = face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
               auto f_phi = face_basis.eval_functions(msh, fc, qp.point());
               auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
               //auto q_f_phi = qp.weight() * f_phi;

               assert(f_phi.size() == fbs);
               assert(c_phi.size() == cell_basis.size());

               for(size_t i = 0; i < fbs; i++){
                  for(size_t j = i; j < fbs; j++){
                     face_mass_matrix(i,j) += qp.weight() * mm_prod(f_phi[i], f_phi[j]);
                  }
               }
               
               //lower part
               for (size_t i = 1; i < fbs; i++)
                  for (size_t j = 0; j < i; j++)
                     face_mass_matrix(i,j) = face_mass_matrix(j,i);

               for(size_t i=0; i< fbs; i++){
                  for(size_t j=0; j< cell_basis.size(); j++){
                     face_trace_matrix(i,j) += qp.weight() * mm_prod(f_phi[i], c_phi[j]);
                  }
               }
            }

            Eigen::LLT<matrix_type> piKF;
            piKF.compute(face_mass_matrix);

            // Step 3a: \pi_F^k( v_F - p_T^k v )
            auto face_range = current_face_range.remove_offset();
            matrix_type MR1 = take(face_trace_matrix, face_range, one_range);
            assert(MR1.cols() == gradrec_oper.rows());
            matrix_type proj2 = piKF.solve(MR1*gradrec_oper);
            matrix_type I_F = matrix_type::Identity(fbs, fbs);
            auto block_offset = current_face_range.min();
            proj2.block(0, block_offset, fbs, fbs) -= I_F;

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


   template<typename Mesh, typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType>
   class elas_like_stabilization_l2
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
      
      elas_like_stabilization_l2()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }
      
      elas_like_stabilization_l2(size_t degree)
      : m_degree(degree)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }
      
      void compute(const mesh_type& msh, const cell_type& cl, const matrix_type& gradrec_oper)
      {
         const size_t DIM= msh.dimension;
         const size_t dpk1 = DIM * binomial(m_degree +1 + DIM, m_degree +1);
         const size_t dpk0 = DIM * binomial( DIM, 0);
         const size_t dpk = DIM * binomial(m_degree  + DIM, m_degree );
         const size_t dpkf = DIM * binomial(m_degree  + DIM -1, m_degree );
         
         assert(cell_basis.size() == dpk1);
          
         auto zero_range         = cell_basis.range(0, m_degree);
         auto one_range          = cell_basis.range(1, m_degree+1);
         
         assert(zero_range.size() == dpk);
         assert(one_range.size() == (dpk1 - dpk0));

         auto fcs = faces(msh, cl);
         auto num_faces = fcs.size();
         
         auto num_cell_dofs = cell_basis.range(0, m_degree).size();
         auto num_face_dofs = face_basis.size();
         
         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);
         
         data = matrix_type::Zero(dsr.total_size(), dsr.total_size());
         
         
         //Step 2: v_T 
         matrix_type I_T = matrix_type::Identity(zero_range.size(), zero_range.size());
         matrix_type proj1  = matrix_type::Zero(zero_range.size(), dsr.total_size());         
         proj1.block(0, 0, zero_range.size(), zero_range.size()) += I_T;
         
         // Step 3: project on faces (eqn. 21)
         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto current_face_range = dsr.face_range(face_i);
            auto fbs = current_face_range.size();
            auto h = diameter(msh, /*fcs[face_i]*/cl);
            auto fc = fcs[face_i];
            
            assert(fbs == dpkf);
            
            matrix_type face_mass_matrix    = matrix_type::Zero(fbs, fbs);
            matrix_type face_trace_matrix   = matrix_type::Zero(fbs, cell_basis.size());
            
            auto face_quadpoints = face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
               auto f_phi = face_basis.eval_functions(msh, fc, qp.point());
               auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
               
               assert(f_phi.size() == fbs);
               assert(c_phi.size() == cell_basis.size());
               
               for(size_t i = 0; i < fbs; i++){
                  for(size_t j = i; j < fbs; j++){
                     face_mass_matrix(i,j) += qp.weight() * mm_prod(f_phi[i], f_phi[j]);
                  }
               }
               
               //lower part
               for (size_t i = 1; i < fbs; i++)
                  for (size_t j = 0; j < i; j++)
                     face_mass_matrix(i,j) = face_mass_matrix(j,i);
                  
                  for(size_t i=0; i< fbs; i++){
                     for(size_t j=0; j< cell_basis.size(); j++){
                        face_trace_matrix(i,j) += qp.weight() * mm_prod(f_phi[i], c_phi[j]);
                     }
                  }
            }
            
            Eigen::LLT<matrix_type> piKF;
            piKF.compute(face_mass_matrix);
            
            // Step 3a: v_F
            auto face_range = current_face_range.remove_offset();

            matrix_type I_F = matrix_type::Identity(fbs, fbs);
            auto block_offset = current_face_range.min();
            matrix_type proj2  = matrix_type::Zero(fbs, dsr.total_size());
            proj2.block(0, block_offset, fbs, fbs) -= I_F;
            
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
   

   template<typename Mesh, //typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType>
   class assembler_elas
   {
      typedef Mesh                                mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;
      typedef typename mesh_type::face            face_type;

      //typedef CellBasisType                       cell_basis_type;
      //typedef CellQuadType                        cell_quadrature_type;
      typedef FaceBasisType                       face_basis_type;
      typedef FaceQuadType                        face_quadrature_type;

      typedef dynamic_matrix<scalar_type>         matrix_type;


      //cell_basis_type                             cell_basis;
      //cell_quadrature_type                        cell_quadrature;

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

      assembler_elas()                 = delete;

      assembler_elas(const mesh_type& msh, size_t degree)
      : m_degree(degree)
      {
         const size_t DIM= msh.dimension;
         const size_t dpkf = DIM * binomial(m_degree  + DIM -1, m_degree );

         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);

         assert(face_basis.size() == dpkf);

         m_num_unknowns = face_basis.size() * (msh.faces_size() + msh.boundary_faces_size());
         matrix = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
         rhs = vector_type::Zero(m_num_unknowns);
      }

      template<typename LocalContrib>
      void
      assemble(const mesh_type& msh, const cell_type& cl, const LocalContrib& lc)
      {

         const size_t DIM= msh.dimension;
         const size_t dpkf = DIM * binomial(m_degree  + DIM -1, m_degree );

         auto fcs = faces(msh, cl);
         std::vector<size_t> l2g(fcs.size() * face_basis.size());

         assert(dpkf == face_basis.size());

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

         assert(lc.first.rows() == fcs.size()*dpkf);
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
      impose_boundary_conditions(const mesh_type& msh, const Function& bc)
      {
         size_t fbs = face_basis.size();
         size_t face_i = 0;
         for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
         {
            auto bfc = *itor;

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first)
               throw std::invalid_argument("This is a bug: face not found");

            size_t face_id = eid.second;

            size_t face_offset = face_id * fbs;
            size_t face_offset_lagrange = (msh.faces_size() + face_i) * fbs;

            auto fqd = face_quadrature.integrate(msh, bfc);

            matrix_type MFF     = matrix_type::Zero(fbs, fbs);
            vector_type rhs_f   = vector_type::Zero(fbs);

            for (auto& qp : fqd)
            {
               auto f_phi = face_basis.eval_functions(msh, bfc, qp.point());

               assert(f_phi.size() == fbs);

               for(size_t i = 0; i < fbs; i++)
                  for(size_t j = i; j < fbs; j++)
                     MFF(i,j) += qp.weight() * mm_prod(f_phi[i], f_phi[j]);
                  
               //lower part
               for (size_t i = 1; i < fbs; i++)
                  for (size_t j = 0; j < i; j++)
                     MFF(i,j) = MFF(j,i);

                  for(size_t i=0; i< fbs; i++)
                     rhs_f(i) += qp.weight() * mm_prod( f_phi[i],  bc(qp.point()));
            }

            #ifdef FILL_COLMAJOR
            for (size_t j = 0; j < MFF.cols(); j++)
            {
               for (size_t i = 0; i < MFF.rows(); i++)
               {
                  m_triplets.push_back( triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i,j)) );
                  m_triplets.push_back( triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i,j)) );
               }
               rhs(face_offset_lagrange+j) = rhs_f(j);
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
   typename FaceBasisType, typename FaceQuadType>
   class projector_elas
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

      projector_elas()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree);
         cell_quadrature     = cell_quadrature_type(2*(m_degree));
         face_basis          = face_basis_type(m_degree+1);
         face_quadrature     = face_quadrature_type(2*(m_degree));
      }

      projector_elas(size_t degree)
      : m_degree(degree)
      {
         cell_basis          = cell_basis_type(m_degree);
         cell_quadrature     = cell_quadrature_type(2*(m_degree));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*(m_degree));
      }

      matrix_type cell_mm;
      matrix_type whole_mm;

      template<typename Function>
      vector_type
      compute_cell(const mesh_type& msh, const cell_type& cl, const Function& f)
      {
         matrix_type mm = matrix_type::Zero(cell_basis.size(), cell_basis.size());
         vector_type rhs = vector_type::Zero(cell_basis.size());

         auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            auto phi = cell_basis.eval_functions(msh, cl, qp.point());

            for(size_t i = 0; i < phi.size(); i++)
               for(size_t j = i; j < phi.size(); j++)
                  mm(i,j)  += qp.weight() *mm_prod(phi[i], phi[j]);
               
            //lower part
            for (size_t i = 1; i < phi.size(); i++)
               for (size_t j = 0; j < i; j++)
                  mm(i,j) = mm(j,i);

               for(size_t i=0; i < phi.size(); i++){
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
         auto fcs = faces(msh, cl);
         vector_type ret = vector_type::Zero(cell_basis.size() + fcs.size()*face_basis.size());
         whole_mm = matrix_type::Zero(cell_basis.size() + fcs.size()*face_basis.size(),
                                      cell_basis.size() + fcs.size()*face_basis.size());

         ret.block(0, 0, cell_basis.size(), 1) = compute_cell(msh, cl, f);
         whole_mm.block(0, 0, cell_basis.size(), cell_basis.size()) = cell_mm;

         size_t face_offset = cell_basis.size();
         for (auto& fc : fcs)
         {
            matrix_type mm = matrix_type::Zero(face_basis.size(), face_basis.size());
            vector_type rhs = vector_type::Zero(face_basis.size());

            auto face_quadpoints = face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
               auto phi = face_basis.eval_functions(msh, fc, qp.point());

               for(size_t i=0; i < phi.size(); i++)
                  for(size_t j=0; j < phi.size(); j++)
                     mm(i,j)  += qp.weight() *mm_prod(phi[i], phi[j]);

                  for(size_t i=0; i < phi.size(); i++)
                     rhs(i) += qp.weight() * mm_prod( f(qp.point()) , phi[i]);

            }

            ret.block(face_offset, 0, face_basis.size(), 1) = mm.llt().solve(rhs);
            whole_mm.block(face_offset, face_offset, face_basis.size(), face_basis.size()) = mm;
            face_offset += face_basis.size();
         }


         return ret;
      }
      
      template<typename Function>
      vector_type
      compute_cell_grad(const mesh_type& msh, const cell_type& cl, const Function& f)
      {
         const size_t DIM = msh.dimension;
         matrix_type mm = matrix_type::Zero(DIM*cell_basis.size(), DIM*cell_basis.size());
         vector_type rhs = vector_type::Zero(DIM*cell_basis.size());
         
         auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            auto phi = cell_basis.eval_functions(msh, cl, qp.point());
            
            for(size_t i = 0; i < phi.size(); i++)
               for(size_t j = i; j < phi.size(); j++)
                  mm(i,j)  += qp.weight() *mm_prod(phi[i], phi[j]);
               
               //lower part
               for (size_t i = 1; i < phi.size(); i++)
                  for (size_t j = 0; j < i; j++)
                     mm(i,j) = mm(j,i);
                  
                  for(size_t i=0; i < phi.size(); i++){
                     rhs(i) += qp.weight() * mm_prod( f(qp.point()) , phi[i]);
                  }
                  
         }
         
         cell_mm = mm;
         return mm.llt().solve(rhs);
      }

   };


} // namespace disk
