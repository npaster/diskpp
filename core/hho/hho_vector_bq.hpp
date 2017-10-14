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
         const size_t face_basis_size = m_bqd.face_basis.range(0, face_degree).size();

         matrix_type stiff_mat = matrix_type::Zero(cell_basis_size, cell_basis_size);

         const auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            const auto dphi = cell_basis.eval_gradients(msh, cl, qp.point());
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
         const auto MG_rowcol_range = cell_basis.range(1, cell_degree + 1);
         const matrix_type MG = take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

         /* RHS, volumetric part. */
         const auto BG_row_range = cell_basis.range(1, cell_degree + 1);
         const auto BG_col_range = cell_basis.range(0, cell_degree);

         const auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();

         const size_t num_cell_dofs = BG_col_range.size();

         const dofspace_ranges dsr(num_cell_dofs, face_basis_size, num_faces);

         assert(dsr.total_size() == (num_cell_dofs + num_faces * face_basis_size));

         matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

         BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
         take(stiff_mat, BG_row_range, BG_col_range);

         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            const auto current_face_range = dsr.face_range(face_i);
            const auto fc = fcs[face_i];
            const auto n = normal(msh, cl, fc);

            const auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);

            for (auto& qp : face_quadpoints)
            {
               auto c_phi = cell_basis.eval_functions(msh, cl, qp.point()); // 0, m_degree);
               const auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point()); // 1, m_degree+1);

               assert(c_phi.size() == cell_basis_size);

               decltype(c_phi) c_dphi_n;

               c_dphi_n.reserve(BG_row_range.to() - BG_row_range.from());

               for(size_t i = BG_row_range.from(); i < BG_row_range.to(); i++){
                  c_dphi_n.push_back(mm_prod(c_dphi[i] , n));
               }

               assert(c_dphi_n.size() == cell_basis_size);
               assert(c_dphi_n.size() == BG.rows());

               matrix_type  T= matrix_type::Zero(BG.rows(), BG_col_range.size());

               for(size_t i = 0; i < BG.rows(); i++){
                  for(size_t j = 0; j < BG_col_range.size(); j++){
                     T(i,j) = qp.weight() * mm_prod(c_dphi_n[i], c_phi[j]);
                  }
               }

               BG.block(0, 0, BG.rows(), BG_col_range.size()) -= T;

               const auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());

               assert(f_phi.size() == face_basis_size);

               matrix_type  F = matrix_type::Zero(BG.rows(), current_face_range.size());

               for(size_t i = 0; i < BG.rows(); i++){
                  for(size_t j = 0; j < current_face_range.size(); j++){
                     F(i,j) = qp.weight() * mm_prod(c_dphi_n[i], f_phi[j]);
                  }
               }

               BG.block(0, current_face_range.min(),
                        BG.rows(), current_face_range.size()) += F;
            }
         }

         assert(MG.rows() == MG.cols());
         assert(MG.cols() == BG.rows());

         oper  = MG.ldlt().solve(BG);    // GT
         data  = BG.transpose() * oper;  // A
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

         template<int DIM>
         std::vector<static_vector<scalar_type, DIM> >
         compute_gphi_n(const std::vector<static_matrix<scalar_type, DIM, DIM> > gphi,
                        const static_vector<scalar_type, DIM> n) const
         {
            std::vector<static_vector<scalar_type, DIM> > ret;

            ret.reserve(gphi.size());

            for(size_t i = 0; i < gphi.size(); i++)
            {
               ret.push_back(mm_prod(gphi[i], n));
            }

            return ret;
         }

         template<int DIM>
         matrix_type
         compute_MG(const std::vector<static_matrix<scalar_type, DIM, DIM> > gphi) const
         {
            const size_t grad_basis_size = gphi.size();
            const size_t grad_degree = m_bqd.grad_degree();
            const size_t poly_space = DIM * DIM * binomial(grad_degree-1 + DIM, grad_degree-1);
            const size_t DIM2 = DIM * DIM;

            matrix_type MG = matrix_type::Zero(grad_basis_size, grad_basis_size);

            // poly classique

            size_t col = 0;
            for(size_t j = 0; j < poly_space; j += DIM2 ){
               for(size_t k = 0; k < DIM; k++ ){//depend de l'ordre des bases
                  for(size_t l = 0; l < DIM; l++ ){//depend de l'ordre des bases
                     for(size_t i = col; i < poly_space; i += DIM2){
                        MG(i,col) = gphi[i](l,k) * gphi[col](l,k);
                     }
                     col++;
                  }
               }
            }

            // RT
            for(std::size_t i = poly_space; i < grad_basis_size; i ++) {
               for(std::size_t j = 0; j < grad_basis_size; j ++) {
                  MG(i,j) = mm_prod(gphi[i], gphi[j]);
               }
            }

            return MG;
         }


         template<int DIM>
         matrix_type
         compute_BG(const std::vector<static_matrix<scalar_type, DIM, DIM> > gphi,
                    const std::vector<static_matrix<scalar_type, DIM, DIM> > cdphi) const
         {
            const size_t grad_basis_size = gphi.size();
            const size_t cell_basis_size = cdphi.size();
            const size_t grad_degree = m_bqd.grad_degree();
            const size_t poly_space = DIM * DIM * binomial(grad_degree-1 + DIM, grad_degree-1);
            const size_t DIM2 = DIM * DIM;

            matrix_type BG = matrix_type::Zero(grad_basis_size, cell_basis_size);

            // poly classique

            size_t row = 0;
            for(size_t i = 0; i < poly_space; i += DIM2 ){
               for(size_t k = 0; k < DIM; k++ ){//depend de l'ordre des bases
                  for(size_t l = 0; l < DIM; l++ ){//depend de l'ordre des bases
                     for(size_t j = l; j < cell_basis_size; j += DIM){
                        BG(row,j) = gphi[row](l,k) * cdphi[j](l,k);
                     }
                     row++;
                  }
               }
            }

            // RT
            for(std::size_t i = poly_space; i < grad_basis_size; i ++) {
               for(std::size_t j = 0; j < cell_basis_size; j ++) {
                  BG(i,j) = mm_prod(gphi[i], cdphi[j]);
               }
            }

            return BG;
         }

      public:

         matrix_type     oper;
         matrix_type     data;

         gradient_reconstruction_elas_full_bq(const BQData& bqd) : m_bqd(bqd)
         {}

         void compute(const mesh_type& msh, const cell_type& cl, const bool compute_data = true)
         {
            const size_t cell_degree = m_bqd.cell_degree();
            const size_t face_degree = m_bqd.face_degree();
            const size_t cell_basis_size = (m_bqd.cell_basis.range(0, cell_degree)).size();
            const size_t face_basis_size = (m_bqd.face_basis.range(0, face_degree)).size();
            const size_t grad_basis_size = m_bqd.grad_basis.size();

            const auto fcs = faces(msh, cl);
            const size_t num_faces = fcs.size();

            timecounter tc;
            double t_base(0.0); double t_cons(0.0); double t_inv(0.0);

            const dofspace_ranges dsr(cell_basis_size, face_basis_size, num_faces);

            assert(dsr.total_size() == (cell_basis_size + num_faces * face_basis_size));

            matrix_type BG = matrix_type::Zero(grad_basis_size, dsr.total_size());

            matrix_type MG = matrix_type::Zero(grad_basis_size, grad_basis_size);


            const auto grad_quadpoints = m_bqd.grad_quadrature.integrate(msh, cl);
            for (auto& qp : grad_quadpoints)
            {
               tc.tic();
               const auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());
               assert(grad_basis_size == gphi.size());

               const auto dphi = m_bqd.cell_basis.eval_gradients(msh, cl, qp.point());
               assert(cell_basis_size == dphi.size());
               tc.toc();
               t_base += tc.to_double();

               tc.tic();

               MG += qp.weight() * compute_MG(gphi);


               BG.block(0,0, grad_basis_size, cell_basis_size) += qp.weight() * compute_BG(gphi, dphi);

               tc.toc();
               t_cons += tc.to_double();
            }

            tc.tic();
            // lower part MG
            for(size_t i = 0; i <  grad_basis_size; i++)
               for(size_t j = i; j < grad_basis_size; j++)
                  MG(i,j) = MG(j,i);

               tc.toc();
            t_cons += tc.to_double();

            for (size_t face_i = 0; face_i < num_faces; face_i++)
            {
               const auto current_face_range = dsr.face_range(face_i);
               const auto fc = fcs[face_i];
               const auto n = normal(msh, cl, fc);
               const auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);

               for (auto& qp : face_quadpoints)
               {
                  tc.tic();
                  const auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point());
                  const auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());
                  tc.toc();
                  t_base += tc.to_double();

                  // tau.n
                  tc.tic();
                  const auto gphi_n = compute_gphi_n(gphi, n);

                  assert(gphi_n.size() == grad_basis_size);

                  matrix_type T = matrix_type::Zero(grad_basis_size, cell_basis_size);

                  for(size_t j = 0; j < cell_basis_size; j++){
                     for(size_t i = 0; i < grad_basis_size; i++){
                        T(i,j) = qp.weight() * mm_prod(gphi_n[i], c_phi[j]);
                     }
                  }

                  assert(T.rows() == grad_basis_size);
                  assert(T.cols() == cell_basis_size);

                  BG.block(0, 0, grad_basis_size, cell_basis_size) -= T;
                  tc.toc();
                  t_cons += tc.to_double();

                  tc.tic();
                  const auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());
                  tc.toc();
                  t_base += tc.to_double();

                  tc.tic();
                  matrix_type F = matrix_type::Zero(BG.rows(), current_face_range.size());

                  for(size_t j = 0; j < current_face_range.size(); j++){
                     for(size_t i = 0; i < grad_basis_size; i++){
                        F(i,j) = qp.weight() * mm_prod(gphi_n[i], f_phi[j]);
                     }
                  }

                  assert(F.rows() == grad_basis_size);
                  assert(F.cols() == current_face_range.size());

                  BG.block(0, current_face_range.min(),
                           grad_basis_size, current_face_range.size()) += F;
                           tc.toc();
                           t_cons += tc.to_double();
               }
            }

            tc.tic();
            oper  = MG.llt().solve(BG);    // GT
            tc.toc();
            t_inv += tc.to_double();
            if(compute_data)
               data  = BG.transpose() * oper;  // A

            assert(oper.rows() == grad_basis_size);
            assert(oper.cols() == dsr.total_size());
         }
      };


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
         const size_t DIM = msh.dimension;
         const size_t cell_degree = m_bqd.cell_degree();
         const size_t face_degree = m_bqd.face_degree();
         const size_t cell_basis_size = (cell_basis.range(0, cell_degree + 1)).size();
         const size_t face_basis_size = (m_bqd.face_basis.range(0, face_degree)).size();

         matrix_type mass_mat = matrix_type::Zero(cell_basis_size, cell_basis_size);

         const auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            const auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
            assert(c_phi.size() == cell_basis_size);

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

         const auto zero_range         = cell_basis.range(0, cell_degree);
         const auto one_range          = cell_basis.range(1, cell_degree+1);

         // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

         //Step 1: compute \pi_T^k p_T^k v (third term).
         const matrix_type M1 = take(mass_mat, zero_range, zero_range);
         const matrix_type M2 = take(mass_mat, zero_range, one_range);
         matrix_type proj1 = -M1.llt().solve(M2*gradrec_oper);

         //Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
         const matrix_type I_T = matrix_type::Identity(zero_range.size(), zero_range.size());
         proj1.block(0, 0, zero_range.size(), zero_range.size()) += I_T;

         const auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();

         const size_t num_cell_dofs = (cell_basis.range(0, cell_degree)).size();

         const dofspace_ranges dsr(num_cell_dofs, face_basis_size, num_faces);

         data = matrix_type::Zero(dsr.total_size(), dsr.total_size());

         // Step 3: project on faces (eqn. 21)
         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            const auto current_face_range = dsr.face_range(face_i);
            const auto h = diameter(msh, /*fcs[face_i]*/cl);
            const auto fc = fcs[face_i];

            matrix_type face_mass_matrix    = matrix_type::Zero(face_basis_size, face_basis_size);
            matrix_type face_trace_matrix   = matrix_type::Zero(face_basis_size, cell_basis_size);

            const auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
               const auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());
               const auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());

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
            const auto face_range = current_face_range.remove_offset();
            const matrix_type MR1 = take(face_trace_matrix, face_range, one_range);
            assert(MR1.cols() == gradrec_oper.rows());

            matrix_type proj2 = piKF.solve(MR1*gradrec_oper);
            const matrix_type I_F = matrix_type::Identity(face_basis_size, face_basis_size);
            const auto block_offset = current_face_range.min();
            proj2.block(0, block_offset, face_basis_size, face_basis_size) -= I_F;

            // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
            const matrix_type MR2 = take(face_trace_matrix, face_range, zero_range);
            assert(MR2.cols() == proj1.rows());
            const matrix_type proj3 = piKF.solve(MR2*proj1);

            assert(proj2.rows() == proj3.rows());
            assert(proj2.cols() == proj3.cols());

            const matrix_type BRF = proj2 + proj3;

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
         const size_t face_basis_size = m_bqd.face_basis.range(0,face_degree).size();

         const auto zero_range = m_bqd.cell_basis.range(0, cell_degree);

         const auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();

         const dofspace_ranges dsr(cell_basis_size, face_basis_size, num_faces);

         data = matrix_type::Zero(dsr.total_size(), dsr.total_size());

         //Step 2: v_T
         const matrix_type I_T = matrix_type::Identity(cell_basis_size, cell_basis_size);
         matrix_type proj1  = matrix_type::Zero(cell_basis_size, dsr.total_size());
         proj1.block(0, 0, cell_basis_size, cell_basis_size) += I_T;

         // Step 3: project on faces (eqn. 21)
         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            const auto current_face_range = dsr.face_range(face_i);
            const auto h = diameter(msh, /*fcs[face_i]*/cl);
            const auto fc = fcs[face_i];

            matrix_type face_mass_matrix    = matrix_type::Zero(face_basis_size, face_basis_size);
            matrix_type face_trace_matrix   = matrix_type::Zero(face_basis_size, cell_basis_size);

            const auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
               const auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());
               const auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point());

               assert(f_phi.size() == face_basis_size);
               assert(c_phi.size() == cell_basis_size);

               for(size_t i = 0; i < face_basis_size; i++){
                  for(size_t j = i; j < face_basis_size; j++){
                     face_mass_matrix(i,j) += qp.weight() * mm_prod(f_phi[i], f_phi[j]);
                  }
               }

               //lower part
               for (size_t i = 1; i < face_basis_size; i++){
                  for (size_t j = 0; j < i; j++){
                     face_mass_matrix(i,j) = face_mass_matrix(j,i);
                  }
               }

               for(size_t i=0; i< face_basis_size; i++){
                  for(size_t j=0; j< cell_basis_size; j++){
                     face_trace_matrix(i,j) += qp.weight() * mm_prod(f_phi[i], c_phi[j]);
                  }
               }
            }

            Eigen::LLT<matrix_type> piKF;
            piKF.compute(face_mass_matrix);

            // Step 3a: v_F
            const auto face_range = current_face_range.remove_offset();

            const matrix_type I_F = matrix_type::Identity(face_basis_size, face_basis_size);
            const auto block_offset = current_face_range.min();
            matrix_type proj2  = matrix_type::Zero(face_basis_size, dsr.total_size());
            proj2.block(0, block_offset, face_basis_size, face_basis_size) -= I_F;

            // Step 3b: \pi_F^k( v_T )
            const matrix_type MR2 = take(face_trace_matrix, face_range, zero_range);
            assert(MR2.cols() == proj1.rows());
            const matrix_type proj3 = piKF.solve(MR2*proj1);

            assert(proj2.rows() == proj3.rows());
            assert(proj2.cols() == proj3.cols());

            // \pi_F^k( v_T ) - v_F
            const matrix_type BRF = proj2 + proj3;

            data += BRF.transpose() * face_mass_matrix * BRF / h;
         }
      }
   };


   template<typename BQData>
      class projector_elas_bq
      {
         typedef typename BQData::mesh_type          mesh_type;
         typedef typename mesh_type::scalar_type     scalar_type;
         typedef typename mesh_type::cell            cell_type;

         typedef typename BQData::cell_basis_type    cell_basis_type;
         typedef typename BQData::cell_quad_type     cell_quad_type;

         typedef dynamic_matrix<scalar_type>         matrix_type;
         typedef dynamic_vector<scalar_type>         vector_type;

         const BQData&                               m_bqd;

      public:

         projector_elas_bq(const BQData& bqd) : m_bqd(bqd)
         {}

         matrix_type cell_mm;
         matrix_type whole_mm;
         matrix_type grad_mm;
         matrix_type pot_mm;

         template<typename Function>
         vector_type
         compute_cell(const mesh_type& msh, const cell_type& cl, const Function& f)
         {
            const size_t cell_degree = m_bqd.cell_degree();
            const size_t cell_basis_size = (m_bqd.cell_basis.range(0, cell_degree)).size();

            matrix_type mm = matrix_type::Zero(cell_basis_size, cell_basis_size);
            vector_type rhs = vector_type::Zero(cell_basis_size);

            const auto cell_quadpoints = m_bqd.cell_quadrature.integrate(msh, cl);
            for (auto& qp : cell_quadpoints)
            {
               const auto phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point());

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
            const size_t face_basis_size = m_bqd.face_basis.range(0,face_degree).size();
            const auto fcs = faces(msh, cl);

            this->compute_cell(msh, cl, f);

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
            const size_t grad_basis_size = m_bqd.grad_basis.size();

            matrix_type mm = matrix_type::Zero(grad_basis_size, grad_basis_size);
            vector_type rhs = vector_type::Zero(grad_basis_size);

            const auto grad_quadpoints = m_bqd.grad_quadrature.integrate(msh, cl);
            for (auto& qp : grad_quadpoints)
            {
               auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());

               for(size_t j = 0; j < grad_basis_size; j++){
                  for(size_t i = j; i < grad_basis_size; i++){
                     mm(i,j) += qp.weight() * mm_prod(gphi[i], gphi[j]);
                  }
               }

               //lower part
               for (size_t i = 0; i < grad_basis_size; i++)
                  for (size_t j = i; j < grad_basis_size; j++)
                     mm(i,j) = mm(j,i);

               for(size_t i=0; i < grad_basis_size; i++){
                  rhs(i) += qp.weight() * mm_prod( f(qp.point()) , gphi[i]);
               }

            }

            grad_mm = mm;
            return mm.llt().solve(rhs);
         }


         template<typename Function>
         vector_type
         compute_pot(const mesh_type& msh, const cell_type& cl, const Function& f)
         {
            const size_t DIM = msh.dimension;
            const size_t cell_degree = m_bqd.cell_degree();
            const size_t cell_basis_size = DIM * binomial(cell_degree +1 +DIM, cell_degree+1);

            cell_basis_type cell_basis        = cell_basis_type(cell_degree + 1);
            cell_quad_type cell_quadrature     = cell_quad_type(2 * (cell_degree + 1));

            assert(cell_basis.size() == cell_basis_size);

            matrix_type mm = matrix_type::Zero(cell_basis_size, cell_basis_size);
            vector_type rhs = vector_type::Zero(cell_basis_size);

            const auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
            for (auto& qp : cell_quadpoints)
            {

               const auto dphi = cell_basis.eval_gradients(msh, cl, qp.point());
               assert(cell_basis_size == dphi.size());

               for(size_t j = 0; j < cell_basis_size; j++){
                  for(size_t i = j; i < cell_basis_size; i++){
                     mm(i,j) += qp.weight() * mm_prod(dphi[i],dphi[j]);
                  }
               }

               for(size_t i=0; i < cell_basis_size; i++){
                  rhs(i) += qp.weight() * mm_prod( f(qp.point()) , dphi[i]);
               }
            }

            // lower part
            for(size_t i = 0; i < cell_basis_size; i++)
               for(size_t j = i; j < cell_basis_size; j++)
                  mm(i,j) = mm(j,i);

            /* LHS: take basis functions derivatives from degree 1 to K+1 */
            const auto MG_rowcol_range = cell_basis.range(1, cell_degree + 1);
            pot_mm = take(mm, MG_rowcol_range, MG_rowcol_range);

            //std::cout << "mm " << mm << '\n';
            //std::cout << "r " << rhs<< '\n';
            return pot_mm.ldlt().solve(rhs.tail(pot_mm.cols()));
         }

      };

      template<typename BQData>
         class displacement_reconstruction_elas_bq
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

            displacement_reconstruction_elas_bq(const BQData& bqd) : m_bqd(bqd)
            {
               cell_basis          = cell_basis_type(bqd.cell_degree() + 1);
               cell_quadrature     = cell_quad_type(2 * (bqd.cell_degree() + 1));
            }


            void compute(const mesh_type& msh, const cell_type& cl)
            {
               const size_t cell_degree = m_bqd.cell_degree();
               const size_t face_degree = m_bqd.face_degree();
               const size_t cell_basis_size = (cell_basis.range(0, cell_degree + 1)).size();
               const size_t face_basis_size = (m_bqd.face_basis.range(0, face_degree)).size();

               matrix_type stiff_mat = matrix_type::Zero(cell_basis_size, cell_basis_size);

               const auto cell_quadpoints = cell_quadrature.integrate(msh, cl);

               for (auto& qp : cell_quadpoints)
               {
                  const auto dphi = cell_basis.eval_gradients(msh, cl, qp.point());
                  assert(cell_basis_size == dphi.size());

                  for(size_t j = 0; j < cell_basis_size; j++){
                     for(size_t i = j; i < cell_basis_size; i++){
                        stiff_mat(i,j) += qp.weight() * mm_prod(dphi[i], dphi[j]);
                     }
                  }
               }

               // lower part
               for(size_t i = 0; i < cell_basis_size; i++)
                  for(size_t j = i; j < cell_basis_size; j++)
                     stiff_mat(i,j) = stiff_mat(j,i);

                  /* LHS: take basis functions derivatives from degree 1 to K+1 */
               const auto MG_rowcol_range = cell_basis.range(1, cell_degree + 1);
               const matrix_type MG = take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

               /* RHS, volumetric part. */
               const auto BG_row_range = cell_basis.range(1, cell_degree + 1);
               const auto BG_col_range = cell_basis.range(0, cell_degree);

               const auto fcs = faces(msh, cl);
               const size_t num_faces = fcs.size();

               const size_t num_cell_dofs = BG_col_range.size();

               const dofspace_ranges dsr(num_cell_dofs, face_basis_size, num_faces);

               assert(dsr.total_size() == (num_cell_dofs + num_faces *face_basis_size));

               matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

               BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
               take(stiff_mat, BG_row_range, BG_col_range);

               for (size_t face_i = 0; face_i < num_faces; face_i++)
               {
                  const auto current_face_range = dsr.face_range(face_i);
                  const auto fc = fcs[face_i];
                  const auto n = normal(msh, cl, fc);

                  const auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);

                  for (auto& qp : face_quadpoints)
                  {
                     auto c_phi = cell_basis.eval_functions(msh, cl, qp.point()); // 0, m_degree);
                     const auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point()); // 1, m_degree+1);

                     decltype(c_phi) c_dphi_n;

                     c_dphi_n.reserve(BG_row_range.to() - BG_row_range.from());

                     for(size_t i=BG_row_range.from(); i< BG_row_range.to(); i++){
                        c_dphi_n.push_back(mm_prod(c_dphi[i] , n));
                     }

                     matrix_type  T= matrix_type::Zero(BG.rows(), BG_col_range.size());

                     assert(c_dphi_n.size() == BG.rows());

                     for(size_t j=0; j<BG_col_range.size(); j++){
                        for(size_t i=0; i< BG.rows(); i++){
                           T(i,j) = qp.weight() * mm_prod(c_dphi_n[i], c_phi[j]);
                        }
                     }

                     BG.block(0, 0, BG.rows(), BG_col_range.size()) -= T;

                     const auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());

                     assert(f_phi.size() == face_basis_size);

                     matrix_type  F = matrix_type::Zero(BG.rows(), current_face_range.size());

                     for(size_t j=0; j < current_face_range.size(); j++){
                        for(size_t i=0; i< BG.rows(); i++){
                           F(i,j) = qp.weight() * mm_prod(c_dphi_n[i],  f_phi[j]);
                        }
                     }

                     BG.block(0, current_face_range.min(),
                              BG.rows(), current_face_range.size()) += F;
                  }
               }

               assert(MG.rows() == MG.cols());
               assert(MG.cols() == BG.rows());

               oper  = MG.ldlt().solve(BG);    // GT
               data  = BG.transpose() * oper;  // A
            }
         };


} // namespace disk
