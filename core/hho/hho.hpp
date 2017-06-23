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


   template<typename CellBasisType, typename CellQuadType, typename Mesh,
   typename Function>
   dynamic_vector<typename Mesh::scalar_type>
   compute_rhs(const Mesh& msh, const typename Mesh::cell& cl,
               const Function& f, size_t degree)
   {
      typedef dynamic_vector<typename Mesh::scalar_type> vector_type;

      auto cell_basis     = CellBasisType(degree);
      auto cell_quad      = CellQuadType(2*degree);

      vector_type ret = vector_type::Zero(cell_basis.size());

      auto cell_quadpoints = cell_quad.integrate(msh, cl);
      for (auto& qp : cell_quadpoints)
      {
         auto phi = cell_basis.eval_functions(msh, cl, qp.point());
         auto fval = f(qp.point());
         for (size_t i = 0; i < cell_basis.size(); i++)
            ret(i) += qp.weight() * mm_prod(fval, phi[i]);
      }

      return ret;
   }

   template<typename Mesh, typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType>
   class projector
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

      projector()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree);
         cell_quadrature     = cell_quadrature_type(2*m_degree);
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }

      projector(size_t degree)
      : m_degree(degree)
      {
         cell_basis          = cell_basis_type(m_degree);
         cell_quadrature     = cell_quadrature_type(2*m_degree);
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }

      matrix_type cell_mm;

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

            mm  += qp.weight() * phi * phi.transpose();
            rhs += qp.weight() * f(qp.point()) * phi;
         }

         cell_mm = mm;
         return mm.llt().solve(rhs);
      }

      template<typename Function>
      vector_type
      compute_whole(const mesh_type& msh, const cell_type& cl, const Function& f)
      {
         auto fcs = faces(msh, cl);
         vector_type ret(cell_basis.size() + fcs.size()*face_basis.size());

         ret.block(0, 0, cell_basis.size(), 1) = compute_cell(msh, cl, f);

         size_t face_offset = cell_basis.size();
         for (auto& fc : fcs)
         {
            matrix_type mm = matrix_type::Zero(face_basis.size(), face_basis.size());
            vector_type rhs = vector_type::Zero(face_basis.size());

            auto face_quadpoints = face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
               auto phi = face_basis.eval_functions(msh, fc, qp.point());

               mm  += qp.weight() * phi * phi.transpose();
               rhs += qp.weight() * f(qp.point()) * phi;
            }

            ret.block(face_offset, 0, face_basis.size(), 1) = mm.llt().solve(rhs);
            face_offset += face_basis.size();
         }

         return ret;
      }
   };

   template<typename BQData>
   class projector_bq
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;
      typedef typename mesh_type::face            face_type;
      typedef typename BQData::cell_basis_type    cell_basis_type;
      typedef typename BQData::face_basis_type    face_basis_type;
      typedef typename BQData::cell_quad_type     cell_quadrature_type;
      typedef typename BQData::face_quad_type     face_quadrature_type;
      typedef dynamic_matrix<scalar_type>         matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;

      const BQData&                               m_bqd;

   public:

      projector_bq(const BQData& bqd) : m_bqd(bqd)
      {}

      matrix_type cell_mm;

      template<typename Function>
      vector_type
      compute_cell(const mesh_type& msh, const cell_type& cl, const Function& f)
      {
         auto cell_degree = m_bqd.cell_degree();
         auto cell_basis_size = m_bqd.cell_basis.range(0, cell_degree).size();

         matrix_type mm = matrix_type::Zero(cell_basis_size, cell_basis_size);
         vector_type rhs = vector_type::Zero(cell_basis_size);

         auto cell_quadpoints = m_bqd.cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            auto phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);

            mm  += qp.weight() * phi * phi.transpose();
            rhs += qp.weight() * f(qp.point()) * phi;
         }

         cell_mm = mm;

         return mm.llt().solve(rhs);
      }

      template<typename Function>
      vector_type
      compute_whole(const mesh_type& msh, const cell_type& cl, const Function& f)
      {
         auto cell_degree = m_bqd.cell_degree();
         auto face_degree = m_bqd.face_degree();
         auto cell_basis_size = m_bqd.cell_basis.range(0, cell_degree).size();
         auto face_basis_size = m_bqd.face_basis.range(0, face_degree).size();


         auto fcs = faces(msh, cl);
         vector_type ret(cell_basis_size + fcs.size()*face_basis_size);

         ret.block(0, 0, cell_basis_size, 1) = compute_cell(msh, cl, f);

         size_t face_offset = cell_basis_size;
         for (auto& fc : fcs)
         {
            matrix_type mm = matrix_type::Zero(face_basis_size, face_basis_size);
            vector_type rhs = vector_type::Zero(face_basis_size);

            auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
               auto phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());

               mm  += qp.weight() * phi * phi.transpose();
               rhs += qp.weight() * f(qp.point()) * phi;
            }

            ret.block(face_offset, 0, face_basis_size, 1) = mm.llt().solve(rhs);
            face_offset += face_basis_size;
         }

         return ret;
      }
   };

   /*
    * struct scalar_potential_gradient_reconstruction;
    *
    * template<typename Mesh, typename What>
    * struct bases_quadratures_trait;
    *
    * template<typename Mesh>
    * struct bases_quadratures_trait<Mesh, scalar_potential_gradient_reconstruction>
    * {
    *    typedef Mesh                                mesh_type;
    *    typedef typename mesh_type::cell            cell_type;
    *    typedef typename mesh_type::face            face_type;
    *    typedef quadrature<mesh_type, cell_type>    cell_quadrature_type;
    *    typedef quadrature<mesh_type, face_type>    face_quadrature_type;
    *    typedef scaled_monomial_scalar_basis<mesh_type, cell_type>  cell_space_basis_type;
    *    typedef scaled_monomial_scalar_basis<mesh_type, face_type>  face_space_basis_type;
    *    typedef scaled_monomial_scalar_basis<mesh_type, cell_type>  result_space_basis_type;
};
*/

   template<typename Mesh,
   template<typename, typename> class Basis,
   template<typename, typename> class Quadrature>
   class basis_quadrature_data /* this name really sucks */
   {
   public:
      typedef Mesh                            mesh_type;
      typedef typename mesh_type::cell        cell_type;
      typedef typename mesh_type::face        face_type;

      typedef Basis<mesh_type, cell_type>         cell_basis_type;
      typedef Basis<mesh_type, face_type>         face_basis_type;
      typedef Quadrature<mesh_type, cell_type>    cell_quad_type;
      typedef Quadrature<mesh_type, face_type>    face_quad_type;

      cell_basis_type     cell_basis;
      face_basis_type     face_basis;
      cell_quad_type      cell_quadrature;
      face_quad_type      face_quadrature;

   private:
      size_t  m_cell_degree, m_face_degree;

      void init(void)
      {
         cell_basis          = cell_basis_type(m_cell_degree+1);
         face_basis          = face_basis_type(m_face_degree);
         cell_quadrature     = cell_quad_type(2*(m_cell_degree+1));
         face_quadrature     = face_quad_type(2*m_face_degree);
      }

   public:
      basis_quadrature_data() : m_cell_degree(1), m_face_degree(1)
      {
         init();
      }

      basis_quadrature_data(size_t cell_degree, size_t face_degree)
      {
         if ( (cell_degree + 1 < face_degree) or (cell_degree > face_degree + 1) )
            throw std::invalid_argument("Invalid cell degree");

         m_cell_degree = cell_degree;
         m_face_degree = face_degree;

         init();
      }

      size_t cell_degree(void) const { return m_cell_degree; }
      size_t face_degree(void) const { return m_face_degree; }
   };

   template<typename BQData>
   class gradient_reconstruction_bq
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;
      typedef typename mesh_type::face            face_type;
      typedef typename BQData::cell_basis_type    cell_basis_type;
      typedef typename BQData::face_basis_type    face_basis_type;
      typedef typename BQData::cell_quad_type     cell_quadrature_type;
      typedef typename BQData::face_quad_type     face_quadrature_type;
      typedef dynamic_matrix<scalar_type>         matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;

      typedef material_tensor<scalar_type, mesh_type::dimension, mesh_type::dimension>
      material_tensor_type;

      const BQData&                               m_bqd;

   public:
      matrix_type     oper;
      matrix_type     data;

      gradient_reconstruction_bq(const BQData& bqd) : m_bqd(bqd)
      {}

      void compute(const mesh_type& msh, const cell_type& cl)
      {
         material_tensor_type id_tens;
         id_tens = material_tensor_type::Identity();
         compute(msh, cl, id_tens);
      }

      void compute(const mesh_type& msh, const cell_type& cl,
                   const material_tensor_type& mtens)
      {
         auto cell_basis_size = m_bqd.cell_basis.size();
         auto face_basis_size = m_bqd.face_basis.size();
         auto cell_degree = m_bqd.cell_degree();

         matrix_type stiff_mat = matrix_type::Zero(cell_basis_size, cell_basis_size);

         auto cell_quadpoints = m_bqd.cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            matrix_type dphi = m_bqd.cell_basis.eval_gradients(msh, cl, qp.point());
            stiff_mat += qp.weight() * dphi * /*mtens **/ dphi.transpose();
         }

         /* LHS: take basis functions derivatives from degree 1 to K+1 */
         auto MG_rowcol_range = m_bqd.cell_basis.range(1, cell_degree+1);
         matrix_type MG = take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

         /* RHS, volumetric part. */
         auto BG_row_range = m_bqd.cell_basis.range(1, cell_degree+1);
         auto BG_col_range = m_bqd.cell_basis.range(0, cell_degree);

         auto fcs = faces(msh, cl);
         auto num_faces = fcs.size();

         auto num_cell_dofs = m_bqd.cell_basis.range(0, cell_degree).size();
         auto num_face_dofs = face_basis_size;

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

         BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
         take(stiff_mat, BG_row_range, BG_col_range);

         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto current_face_range = dsr.face_range(face_i);
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);
            auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);

            for (auto& qp : face_quadpoints)
            {
               matrix_type c_phi =
               m_bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
               matrix_type c_dphi =
               m_bqd.cell_basis.eval_gradients(msh, cl, qp.point(), 1, cell_degree+1);

               matrix_type c_dphi_n = (c_dphi /** mtens*/ * n);
               matrix_type T = qp.weight() * c_dphi_n * c_phi.transpose();

               BG.block(0, 0, BG.rows(), BG_col_range.size()) -= T;

               matrix_type f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());
               matrix_type F = qp.weight() * c_dphi_n * f_phi.transpose();

               BG.block(0, current_face_range.min(),
                        BG.rows(), current_face_range.size()) += F;
            }
         }

         oper  = MG.llt().solve(BG);    // GT
         data  = BG.transpose() * oper;  // A
      }

   };

   template<typename Mesh, typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType>
   class gradient_reconstruction
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

      gradient_reconstruction()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }

      gradient_reconstruction(size_t degree)
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
         matrix_type stiff_mat = matrix_type::Zero(cell_basis.size(), cell_basis.size());

         auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            matrix_type dphi = cell_basis.eval_gradients(msh, cl, qp.point());
            stiff_mat += qp.weight() * dphi * (/*mtens **/ dphi.transpose());
         }

         /* LHS: take basis functions derivatives from degree 1 to K+1 */
         auto MG_rowcol_range = cell_basis.range(1, m_degree+1);
         matrix_type MG = take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

         /* RHS, volumetric part. */
         auto BG_row_range = cell_basis.range(1, m_degree+1);
         auto BG_col_range = cell_basis.range(0, m_degree);

         auto fcs = faces(msh, cl);
         auto num_faces = fcs.size();

         auto num_cell_dofs = cell_basis.range(0, m_degree).size();
         auto num_face_dofs = face_basis.size();

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

         BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
         take(stiff_mat, BG_row_range, BG_col_range);

         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto current_face_range = dsr.face_range(face_i);
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);
            auto face_quadpoints = face_quadrature.integrate(msh, fc);

            for (auto& qp : face_quadpoints)
            {
               matrix_type c_phi =
               cell_basis.eval_functions(msh, cl, qp.point(), 0, m_degree);
               matrix_type c_dphi =
               cell_basis.eval_gradients(msh, cl, qp.point(), 1, m_degree+1);

               matrix_type c_dphi_n = (c_dphi /** mtens*/) * n;
               matrix_type T = qp.weight() * c_dphi_n * c_phi.transpose();

               BG.block(0, 0, BG.rows(), BG_col_range.size()) -= T;

               matrix_type f_phi = face_basis.eval_functions(msh, fc, qp.point());
               matrix_type F = qp.weight() * c_dphi_n * f_phi.transpose();

               BG.block(0, current_face_range.min(),
                        BG.rows(), current_face_range.size()) += F;
            }
         }

         oper  = MG.llt().solve(BG);    // GT
         data  = BG.transpose() * oper;  // A
      }
   };



   template<typename Mesh, typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType,
   typename GradBasisType, typename GradQuadType>
   class gradient_reconstruction_full
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

      typedef material_tensor<scalar_type, mesh_type::dimension, mesh_type::dimension>
      material_tensor_type;

      cell_basis_type                             cell_basis;
      cell_quadrature_type                        cell_quadrature;

      face_basis_type                             face_basis;
      face_quadrature_type                        face_quadrature;

      grad_basis_type                             grad_basis;
      grad_quadrature_type                        grad_quadrature;

      size_t                                      m_degree;

   public:
      matrix_type     oper;
      matrix_type     data;

      gradient_reconstruction_full()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree);
         cell_quadrature     = cell_quadrature_type(2*m_degree);
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
         grad_basis          = grad_basis_type(m_degree);
         grad_quadrature     = grad_quadrature_type(2*m_degree);
      }

      gradient_reconstruction_full(size_t degree)
      : m_degree(degree)
      {
         cell_basis          = cell_basis_type(m_degree);
         cell_quadrature     = cell_quadrature_type(2*m_degree);
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
         grad_basis          = grad_basis_type(m_degree);
         grad_quadrature     = grad_quadrature_type(2*m_degree);
      }

      void compute(const mesh_type& msh, const cell_type& cl)
      {
         auto BG_col_range = cell_basis.range(0, m_degree);

         const size_t  DIM = msh.dimension;
         const size_t gpk = DIM * BG_col_range.size();

         auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();

         const size_t num_grad_dofs = grad_basis.range(0, m_degree).size();
         const size_t num_cell_dofs = cell_basis.range(0, m_degree).size();
         const size_t num_face_dofs = face_basis.size();

         assert(num_grad_dofs == gpk);

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         assert(dsr.total_size() == (num_cell_dofs + num_faces * num_face_dofs));

         matrix_type BG = matrix_type::Zero(num_grad_dofs, dsr.total_size());

         matrix_type MG = matrix_type::Zero(num_grad_dofs, num_grad_dofs);


         auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            auto gphi = grad_basis.eval_functions(msh, cl, qp.point());
            assert(num_grad_dofs == gphi.size());

            auto dphi = cell_basis.eval_gradients(msh, cl, qp.point(), 0, m_degree);
            assert(cell_basis.size() == dphi.rows());

            for(size_t i = 0; i < gphi.size(); i++){
               for(size_t j = i; j < gphi.size(); j++){
                  MG(i,j) += qp.weight() * mm_prod(gphi[i], gphi[j]);
               }

               for(size_t j = BG_col_range.from(); j < BG_col_range.to(); j++){
                  auto dphi_j = gphi[i] ;
                  for(size_t d = 0; d < dphi.cols(); d++) { dphi_j(d) = dphi(j,d);}
                     BG(i,j) += qp.weight() * mm_prod(gphi[i], dphi_j );
               }
            }
         }


         // lower part MG
         for(size_t i = 1; i < gpk; i++)
            for(size_t j = 0; j < i; j++)
               MG(i,j) = MG(j,i);

         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto current_face_range = dsr.face_range(face_i);
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);
            auto face_quadpoints = face_quadrature.integrate(msh, fc);

            for (auto& qp : face_quadpoints)
            {
               matrix_type c_phi = cell_basis.eval_functions(msh, cl, qp.point(), 0, m_degree);
               auto gphi = grad_basis.eval_functions(msh, cl, qp.point());

               // tau.n
               matrix_type gphi_n = matrix_type::Zero(gphi.size(), 1);
               for(size_t i = 0; i < gphi.size(); i++)
                  gphi_n(i,0) = mm_prod(gphi[i], n);


               matrix_type T = qp.weight() * gphi_n * c_phi.transpose();

               assert(T.rows() == BG.rows());
               assert(T.cols() == BG_col_range.size());

               BG.block(0, 0, BG.rows(), BG_col_range.size()) -= T;

               matrix_type f_phi = face_basis.eval_functions(msh, fc, qp.point());
               matrix_type F = qp.weight() * gphi_n * f_phi.transpose();

               assert(F.rows() == BG.rows());
               assert(F.cols() == current_face_range.size());

               BG.block(0, current_face_range.min(),
                        BG.rows(), current_face_range.size()) += F;
            }
         }

         oper  = MG.llt().solve(BG);    // GT
         data  = BG.transpose() * oper;  // A
      }
   };

   template<typename Mesh, typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType,
   typename DivCellBasisType, typename DivCellQuadType>
   class divergence_reconstruction
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

      divergence_reconstruction()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
         div_cell_basis      = div_cell_basis_type(m_degree);
         div_cell_quadrature = div_cell_quadrature_type(2*m_degree);
      }

      divergence_reconstruction(size_t degree)
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
         auto dcbs = div_cell_basis.size();
         matrix_type MD = matrix_type::Zero(dcbs, dcbs);

         auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            auto phi_d = div_cell_basis.eval_functions(msh, cl, qp.point());
            for (size_t i = 0; i < dcbs; i++)
               for (size_t j = 0; j < dcbs; j++)
                  MD(i,j) += qp.weight() * mm_prod(phi_d[i], phi_d[j]);
         }

         /* RHS, volumetric part. */
         auto fcs = faces(msh, cl);
         auto num_cell_dofs = cell_basis.range(0, m_degree).size();
         auto num_face_dofs = face_basis.size();
         auto num_faces = fcs.size();

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         matrix_type BD = matrix_type::Zero(dcbs, dsr.total_size());
         for (auto& qp : cell_quadpoints)
         {
            auto phi = cell_basis.eval_functions(msh, cl, qp.point());
            auto dphi_d = div_cell_basis.eval_gradients(msh, cl, qp.point());

            for (size_t i = 0; i < dphi_d.size(); i++)
               for (size_t j = 0; j < num_cell_dofs; j++)
                  BD(i,j) -= qp.weight() * mm_prod(dphi_d[i], phi[j]);
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
               for (size_t i = 0; i < phi_d.size(); i++)
               {
                  for (size_t j = 0; j < face_basis.size(); j++)
                  {
                     auto p1 = mm_prod(phi[j], n);
                     scalar_type p2 = mm_prod(p1, phi_d[i]);
                     BD(i,face_offset+j) += qp.weight() * p2;
                  }
               }
            }

            face_offset += face_basis.size();
         }

         oper = MD.partialPivLu().solve(BD);
         data = BD.transpose() * oper;
      }
   };




   template<typename BQData>
   class diffusion_like_stabilization_bq
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;
      typedef typename mesh_type::face            face_type;
      typedef typename BQData::cell_basis_type    cell_basis_type;
      typedef typename BQData::face_basis_type    face_basis_type;
      typedef typename BQData::cell_quad_type     cell_quadrature_type;
      typedef typename BQData::face_quad_type     face_quadrature_type;
      typedef dynamic_matrix<scalar_type>         matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;

      const BQData&                               m_bqd;

   public:
      matrix_type     data;

      diffusion_like_stabilization_bq(const BQData& bqd) : m_bqd(bqd)
      {}

      void compute(const mesh_type& msh, const cell_type& cl, const matrix_type& gradrec_oper)
      {
         auto cell_basis_size = m_bqd.cell_basis.size();
         auto face_basis_size = m_bqd.face_basis.size();
         auto cell_degree = m_bqd.cell_degree();

         matrix_type mass_mat = matrix_type::Zero(cell_basis_size, cell_basis_size);

         auto cell_quadpoints = m_bqd.cell_quadrature.integrate(msh, cl);
         for (auto& qp : cell_quadpoints)
         {
            auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point());
            mass_mat += qp.weight() * c_phi * c_phi.transpose();
         }

         auto zero_range         = m_bqd.cell_basis.range(0, cell_degree);
         auto one_range          = m_bqd.cell_basis.range(1, cell_degree+1);

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

         auto num_cell_dofs = m_bqd.cell_basis.range(0, cell_degree).size();
         auto num_face_dofs = face_basis_size;

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         data = matrix_type::Zero(dsr.total_size(), dsr.total_size());

         // Step 3: project on faces (eqn. 21)
         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto current_face_range = dsr.face_range(face_i);
            auto fbs = current_face_range.size();
            auto h = diameter(msh, /*fcs[face_i]*/cl);
            auto fc = fcs[face_i];

            matrix_type face_mass_matrix    = matrix_type::Zero(fbs, fbs);
            matrix_type face_trace_matrix   = matrix_type::Zero(fbs, cell_basis_size);

            auto face_quadpoints = m_bqd.face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
               auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());
               auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point());
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
   typename FaceBasisType, typename FaceQuadType>
   class diffusion_like_stabilization
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

      diffusion_like_stabilization()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }

      diffusion_like_stabilization(size_t degree)
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

   template<typename BQData>
   class diffusion_like_static_condensation_bq
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;
      typedef typename mesh_type::face            face_type;
      typedef typename BQData::cell_basis_type    cell_basis_type;
      typedef typename BQData::face_basis_type    face_basis_type;
      typedef typename BQData::cell_quad_type     cell_quadrature_type;
      typedef typename BQData::face_quad_type     face_quadrature_type;
      typedef dynamic_matrix<scalar_type>         matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;

      const BQData&                               m_bqd;

   public:
      diffusion_like_static_condensation_bq(const BQData& bqd) : m_bqd(bqd)
      {}

      auto
      compute(const mesh_type& msh, const cell_type& cl,
              const matrix_type& local_mat,
              const vector_type& cell_rhs)
      {
         auto cell_degree = m_bqd.cell_degree();
         auto face_degree = m_bqd.face_degree();
         auto num_cell_dofs = m_bqd.cell_basis.range(0, cell_degree).size();
         auto num_face_dofs = m_bqd.face_basis.range(0, face_degree).size();

         auto fcs = faces(msh, cl);
         auto num_faces = fcs.size();

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         size_t cell_size = dsr.cell_range().size();
         size_t face_size = dsr.all_faces_range().size();

         matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
         matrix_type K_TF = local_mat.topRightCorner(cell_size, face_size);
         matrix_type K_FT = local_mat.bottomLeftCorner(face_size, cell_size);
         matrix_type K_FF = local_mat.bottomRightCorner(face_size, face_size);

         assert(K_TT.cols() == cell_size);
         assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
         assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
         assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
         assert(K_FT.cols() + K_FF.cols() == local_mat.cols());

         auto K_TT_ldlt = K_TT.llt();
         matrix_type AL = K_TT_ldlt.solve(K_TF);
         vector_type bL = K_TT_ldlt.solve(cell_rhs);

         matrix_type AC = K_FF - K_FT * AL;
         vector_type bC = /* no projection on faces, eqn. 26*/ - K_FT * bL;

         return std::make_pair(AC, bC);
      }

      vector_type
      recover(const mesh_type& msh, const cell_type& cl,
              const matrix_type& local_mat,
              const vector_type& cell_rhs,
              const vector_type& solF)
      {
         auto cell_degree = m_bqd.cell_degree();
         auto face_degree = m_bqd.face_degree();
         auto num_cell_dofs = m_bqd.cell_basis.range(0, cell_degree).size();
         auto num_face_dofs = m_bqd.face_basis.range(0, face_degree).size();

         auto fcs = faces(msh, cl);
         auto num_faces = fcs.size();

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         size_t cell_size        = dsr.cell_range().size();
         size_t all_faces_size   = dsr.all_faces_range().size();

         vector_type ret( dsr.total_size() );

         matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
         matrix_type K_TF = local_mat.topRightCorner(cell_size, all_faces_size);

         vector_type solT = K_TT.llt().solve(cell_rhs - K_TF*solF);

         ret.head(cell_size)         = solT;
         ret.tail(all_faces_size)    = solF;

         return ret;
      }

   };


   template<typename Mesh, typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType>
   class diffusion_like_static_condensation
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
      diffusion_like_static_condensation()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree);
         cell_quadrature     = cell_quadrature_type(2*m_degree);
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }

      diffusion_like_static_condensation(size_t degree)
      : m_degree(degree)
      {
         cell_basis          = cell_basis_type(m_degree);
         cell_quadrature     = cell_quadrature_type(2*m_degree);
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }

      std::pair<matrix_type, vector_type>
      compute(const mesh_type& msh, const cell_type& cl,
              const matrix_type& local_mat,
              const vector_type& cell_rhs)
      {
         auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();

         const size_t num_cell_dofs = cell_basis.size();
         const size_t num_face_dofs = face_basis.size();

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         const size_t cell_size = dsr.cell_range().size();
         const size_t face_size = dsr.all_faces_range().size();

         assert( cell_size ==  cell_rhs.rows() && "wrong rhs dimension");
         assert( (cell_size + face_size) == local_mat.rows() && "wrong lhs rows dimension");
         assert( (cell_size + face_size) == local_mat.cols() && "wrong lhs cols dimension");

         matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
         matrix_type K_TF = local_mat.topRightCorner(cell_size, face_size);
         matrix_type K_FT = local_mat.bottomLeftCorner(face_size, cell_size);
         matrix_type K_FF = local_mat.bottomRightCorner(face_size, face_size);

         assert(K_TT.cols() == cell_size);
         assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
         assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
         assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
         assert(K_FT.cols() + K_FF.cols() == local_mat.cols());

         auto K_TT_ldlt = K_TT.llt();
         matrix_type AL = K_TT_ldlt.solve(K_TF);
         vector_type bL = K_TT_ldlt.solve(cell_rhs);

         matrix_type AC = K_FF - K_FT * AL;
         vector_type bC = /* no projection on faces, eqn. 26*/ - K_FT * bL;

         return std::make_pair(AC, bC);
      }

      std::pair<matrix_type, vector_type>
      compute(const mesh_type& msh, const cell_type& cl,
              const matrix_type& local_mat,
              const vector_type& rhs, bool l_face)
      {
         auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();

         const size_t num_cell_dofs = cell_basis.size();
         const size_t num_face_dofs = face_basis.size();

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         const size_t cell_size = dsr.cell_range().size();
         const size_t faces_size = dsr.all_faces_range().size();

         assert( (cell_size + faces_size) ==  rhs.rows() && "wrong rhs dimension");
         assert( (cell_size + faces_size) ==  local_mat.rows() && "wrong lhs rows dimension");
         assert( (cell_size + faces_size) ==  local_mat.cols() && "wrong lhs cols dimension");

         matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
         matrix_type K_TF = local_mat.topRightCorner(cell_size, faces_size);
         matrix_type K_FT = local_mat.bottomLeftCorner(faces_size, cell_size);
         matrix_type K_FF = local_mat.bottomRightCorner(faces_size, faces_size);

         assert(K_TT.cols() == cell_size && "wrong K_TT dimension");
         assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
         assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
         assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
         assert(K_FT.cols() + K_FF.cols() == local_mat.cols());

         vector_type rhs_cell = rhs.block(0,0, cell_size, 1);
         vector_type rhs_faces = rhs.block(cell_size,0, faces_size, 1);

         assert((rhs_cell.rows() + rhs_faces.rows() ) ==  rhs.rows() && "wrong rhs decomposition");

         auto K_TT_ldlt = K_TT.llt();
         matrix_type AL = K_TT_ldlt.solve(K_TF);
         vector_type bL = K_TT_ldlt.solve(rhs_cell);

         matrix_type AC = K_FF - K_FT * AL;
         vector_type bC = rhs_faces - K_FT * bL;

         return std::make_pair(AC, bC);
      }

      vector_type
      recover(const mesh_type& msh, const cell_type& cl,
              const matrix_type& local_mat,
              const vector_type& cell_rhs,
              const vector_type& solF)
      {
         auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();

         const size_t num_cell_dofs = cell_basis.size();
         size_t num_face_dofs = face_basis.size();

         dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

         const size_t cell_size        = dsr.cell_range().size();
         const size_t all_faces_size   = dsr.all_faces_range().size();

         assert(cell_rhs.rows() == cell_size && "wrong rhs_cell dimension");
         assert(solF.rows() == all_faces_size && "wrong solF dimension");
         assert( (cell_size + all_faces_size) == local_mat.rows() && "wrong lhs rows dimension");
         assert( (cell_size + all_faces_size) == local_mat.cols() && "wrong lhs cols dimension");

         vector_type ret( dsr.total_size() );

         matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
         matrix_type K_TF = local_mat.topRightCorner(cell_size, all_faces_size);

         assert(K_TT.cols() == cell_size && "wrong K_TT dimension");
         assert(K_TT.cols() + K_TF.cols() == local_mat.cols());

         assert(num_cell_dofs == cell_rhs.size());
         assert((num_face_dofs * num_faces) == solF.size());

         vector_type solT = K_TT.llt().solve(cell_rhs - K_TF*solF);

         ret.head(cell_size)         = solT;
         ret.tail(all_faces_size)    = solF;

         return ret;
      }



   };

   template<typename Mesh>
   class multiscale_local_problem
   {
      typedef Mesh                                mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;
      typedef typename mesh_type::face            face_type;
      typedef dynamic_matrix<scalar_type>         matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;
      typedef Eigen::SparseMatrix<scalar_type>    sparse_matrix_type;
      typedef Eigen::Triplet<scalar_type>         triplet_type;

      typedef basis_quadrature_data<mesh_type,
      scaled_monomial_scalar_basis,
      quadrature>   bqdata_type;

      size_t                        m_degree,         /* degree requested by the user */
      m_inner_degree,   /* fine scale method degree */
      m_cf_degree,      /* cell function degree (\Phi_T) */
      m_ff_degree;      /* face function degree (\Phi_F) */

      size_t                        funcnum, ptlevels, reflevels;

   public:
      sparse_matrix_type                          matrix;
      matrix_type                                 rhs;

      multiscale_local_problem()
      {
         m_degree = 1;
         m_inner_degree = 2;
         m_cf_degree = 0;
         m_ff_degree = 1;
         load_params();
      }

      multiscale_local_problem(size_t degree)
      {
         // TODO: add checks
         m_degree = degree;
         m_inner_degree = degree+1;
         m_cf_degree = degree-1;
         m_ff_degree = degree;
         load_params();
      }

      void load_params(void)
      {
         // commented for the moment
         int i=0;
         //sol::state lua;
         //lua.script_file("params.lua");
         //funcnum = lua.get<size_t>("funcnum");
         //ptlevels = lua.get<size_t>("ptlevels");
         //reflevels = lua.get<size_t>("reflevels");
      }

      //template<typename OuterCellBasis, OuterFaceBasis>
      void assemble(const mesh_type& coarse_msh,
                    const cell_type& coarse_cl/*,
                    const OuterCellBasis& outer_cell_basis,
                    const OuterFaceBasis& outer_face_basis*/)
      {
         assert(m_degree > 0);
         submesher<Mesh>                                 submesher;
         bqdata_type                                     bqd(m_inner_degree, m_inner_degree);
         gradient_reconstruction_bq<bqdata_type>         gradrec(bqd);
         diffusion_like_stabilization_bq<bqdata_type>    stab(bqd);

         scaled_monomial_scalar_basis<mesh_type, cell_type> outer_cell_basis(m_cf_degree);
         scaled_monomial_scalar_basis<mesh_type, face_type> outer_face_basis(m_ff_degree);

         auto msh = submesher.generate_mesh(coarse_msh, coarse_cl, reflevels);

         auto num_cell_dofs = howmany_dofs(bqd.cell_basis, 0, m_inner_degree);
         auto num_face_dofs = howmany_dofs(bqd.face_basis, 0, m_inner_degree);

         auto num_cell_funcs = howmany_dofs(outer_cell_basis, 0, m_cf_degree);
         auto num_face_funcs = howmany_dofs(outer_face_basis, 0, m_ff_degree);

         auto matrix_face_offset = num_cell_dofs * msh.cells_size();
         auto matrix_mult_offset = matrix_face_offset + num_face_dofs * msh.faces_size();
         auto system_size = num_cell_dofs * msh.cells_size() +
         num_face_dofs * msh.faces_size() +
         num_face_funcs * 3;

         matrix = sparse_matrix_type(system_size, system_size);
         rhs = matrix_type::Zero(system_size, num_cell_funcs);

         std::vector<triplet_type> triplets;

         /* Assemble standard HHO part */
         size_t cell_idx = 0;
         for (auto& cl : msh)
         {
            auto fcs = faces(msh, cl);
            std::vector<size_t> l2g(num_cell_dofs + fcs.size() * num_face_dofs);

            /* Build DOF offset table: cell */
            for (size_t i = 0; i < num_cell_dofs; i++)
               l2g[i] = cell_idx * num_cell_dofs + i;

            /* Build DOF offset table: faces */
            for (size_t i = 0; i < fcs.size(); i++)
            {
               auto fc = fcs[i];
               //std::cout << msh.boundary_id(fc) << std::endl;
               auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
               if (!eid.first)
                  throw std::invalid_argument("This is a bug: face not found");

               auto face_id = eid.second;
               /* global offset of current face */
               auto face_offset = matrix_face_offset + face_id * num_face_dofs;

               /* offset in the DOF table */
               auto dt_ofs = num_cell_dofs + i * num_face_dofs;

               for (size_t j = 0; j < num_face_dofs; j++)
                  l2g[dt_ofs+j] = face_offset+j;
            }

            /* Compute HHO element contribution */
            gradrec.compute(msh, cl);
            stab.compute(msh, cl, gradrec.oper);
            dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
            assert(loc.rows() == l2g.size());
            assert(loc.cols() == l2g.size());

            /* Assemble into the matrix */
            for (size_t i = 0; i < l2g.size(); i++)
               for (size_t j = 0; j < l2g.size(); j++)
                  triplets.push_back( triplet_type(l2g[i], l2g[j], loc(i,j)) );

               /* Now compute the multiscale-specific stuff */
               matrix_type cell_rhs;
            cell_rhs = matrix_type::Zero(num_cell_dofs, num_cell_funcs);

            auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
            for (auto& qp : cell_quadpoints)
            {
               auto phi = bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, m_inner_degree);
               auto phi_coarse = outer_cell_basis.eval_functions(coarse_msh, coarse_cl, qp.point(), 0, m_cf_degree);
               cell_rhs += qp.weight() * phi * phi_coarse.transpose();
            }

            for (size_t i = 0; i < num_cell_dofs; i++)
            {
               auto cell_offset = cell_idx * num_cell_dofs;
               for (size_t j = 0; j < num_cell_funcs; j++)
                  rhs(cell_offset+i, j) = cell_rhs(i,j);
            }

            cell_idx++;
         }

         auto bid_list = msh.boundary_id_list();

         for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
         {
            auto fc = *itor;
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
               throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            /* row */
            auto face_offset = matrix_face_offset + face_id * num_face_dofs;

            size_t coarse_id = msh.boundary_id(fc);
            size_t fine_id;
            for (size_t i = 0; i < bid_list.size(); i++)
               if (bid_list[i] == coarse_id)
               {
                  fine_id = i;
                  break;
               }

               assert (fine_id < 3);

            auto coarse_fc = *(coarse_msh.faces_begin() + coarse_id);

            matrix_type face_mass_matrix;
            face_mass_matrix = matrix_type::Zero(num_face_dofs, num_face_funcs);

            auto face_quadpoints = bqd.face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
               auto phi = bqd.face_basis.eval_functions(msh, fc, qp.point());
               auto phi_coarse = outer_face_basis.eval_functions(msh, coarse_fc, qp.point(), 0, m_ff_degree);
               face_mass_matrix += qp.weight() * phi * phi_coarse.transpose();
            }

            /* col */
            auto mult_offset = matrix_mult_offset + fine_id * num_face_funcs;

            for (size_t j = 0; j < num_face_dofs; j++)
            {
               for (size_t k = 0; k < num_face_funcs; k++)
               {
                  size_t row = face_offset + j;
                  size_t col = mult_offset + k;
                  triplets.push_back( triplet_type(row, col, -face_mass_matrix(j,k)) );
                  triplets.push_back( triplet_type(col, row, -face_mass_matrix(j,k)) );
               }
            }
         }

         matrix.setFromTriplets(triplets.begin(), triplets.end());

         std::stringstream ssw;
         ssw << "matrix.mat";
         std::ofstream ofs(ssw.str());

         for (int k=0; k<matrix.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(matrix,k); it; ++it)
               ofs << it.row() << " " << it.col() << " " << it.value() << std::endl;

            ofs.close();

         triplets.clear();



         #ifdef HAVE_INTEL_MKL
         Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
         #else
         Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
         #endif

         solver.analyzePattern(matrix);
         solver.factorize(matrix);
         matrix_type X = solver.solve(rhs);

         auto eid = find_element_id(coarse_msh.cells_begin(), coarse_msh.cells_end(), coarse_cl);
         if (!eid.first)
            throw std::invalid_argument("This is a bug: cell not found");

         auto cell_id = eid.second;

         std::stringstream ss;
         ss << "multiscale_plot_" << cell_id << ".dat";

         std::ofstream ofs_sol(ss.str());
         cell_idx = 0;
         for (auto& cl : msh)
         {
            auto cell_sol = X.block(num_cell_dofs * cell_idx, funcnum, num_cell_dofs, 1);
            //auto qps = bqd.cell_quadrature.integrate(msh, cl);
            auto qps = make_test_points(msh, cl, ptlevels);
            for (auto& qp : qps)
            {
               auto phi = bqd.cell_basis.eval_functions(msh, cl, qp);

               scalar_type pot = 0.0;
               for (size_t i = 0; i < bqd.cell_basis.range(0, m_inner_degree).size(); i++)
                  pot += phi[i] * cell_sol(i,0);

               auto tp = qp;
               for (size_t i = 0; i < mesh_type::dimension; i++)
                  ofs_sol << tp[i] << " ";
               ofs_sol << pot << std::endl;
            }

            cell_idx++;
         }
         ofs_sol.close();

      }
   };

   template<typename Mesh, //typename CellBasisType, typename CellQuadType,
   typename FaceBasisType, typename FaceQuadType>
   class assembler
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

      assembler()                 = delete;

      assembler(const mesh_type& msh, size_t degree)
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
   class diffusion_like_stabilization_l2
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

      diffusion_like_stabilization_l2()
      : m_degree(1)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }

      diffusion_like_stabilization_l2(size_t degree)
      : m_degree(degree)
      {
         cell_basis          = cell_basis_type(m_degree+1);
         cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
         face_basis          = face_basis_type(m_degree);
         face_quadrature     = face_quadrature_type(2*m_degree);
      }

      void compute(const mesh_type& msh, const cell_type& cl, const matrix_type& gradrec_oper)
      {
         auto zero_range         = cell_basis.range(0, m_degree);
         auto one_range          = cell_basis.range(1, m_degree+1);


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

} // namespace disk
