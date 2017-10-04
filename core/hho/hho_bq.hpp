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
         face_quadrature     = face_quad_type(2 * m_face_degree);
         grad_quadrature     = grad_quad_type(2 * m_grad_degree + 2);
      }

   public:
      basis_quadrature_data_full() : m_cell_degree(1), m_face_degree(1), m_grad_degree(1)
      {
         init();
      }

      basis_quadrature_data_full(const size_t face_degree, const size_t cell_degree, const size_t grad_degree)
      {
         if ( (cell_degree + 1 < face_degree) or (cell_degree > face_degree + 1) )
            throw std::invalid_argument("Invalid cell degree");

         if (grad_degree < cell_degree)
            throw std::invalid_argument("Invalid grad degree");

         m_cell_degree = cell_degree;
         m_face_degree = face_degree;
         m_grad_degree = grad_degree;

         init();
      }

      void info_degree() const
      {
         std::cout << "Face degree: "  << m_face_degree  << std::endl;
         std::cout << "Cell degree: "  << m_cell_degree  << std::endl;
         std::cout << "Grad degree: "  << m_grad_degree  << std::endl;

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

      typedef typename BQData::face_quad_type        face_quadrature_type;

      typedef dynamic_matrix<scalar_type>         matrix_type;
      typedef dynamic_vector<scalar_type>         vector_type;

      const BQData&                               m_bqd;

      face_quadrature_type face_quadrature;

   public:

      matrix_type     oper;
      matrix_type     data;

      gradient_reconstruction_full_bq(const BQData& bqd) : m_bqd(bqd)
      {
         face_quadrature  = face_quadrature_type(m_bqd.face_degree() + m_bqd.grad_degree());
      }

      void compute(const mesh_type& msh, const cell_type& cl, const bool compute_data = true)
      {
         const size_t DIM = msh.dimension;
         const size_t cell_degree = m_bqd.cell_degree();
         const size_t face_degree = m_bqd.face_degree();
         const size_t grad_degree = m_bqd.grad_degree();
         const size_t cell_basis_size = (m_bqd.cell_basis.range(0, cell_degree)).size();
         const size_t face_basis_size = (m_bqd.face_basis.range(0, face_degree)).size();
         const size_t grad_basis_size = (m_bqd.grad_basis.range(0, grad_degree)).size();


         auto fcs = faces(msh, cl);
         const size_t num_faces = fcs.size();

         assert(grad_basis_size == DIM * binomial( grad_degree +DIM, grad_degree));

         dofspace_ranges dsr(cell_basis_size, face_basis_size, num_faces);

         assert(dsr.total_size() == (cell_basis_size + num_faces * face_basis_size));

         matrix_type BG = matrix_type::Zero(grad_basis_size, dsr.total_size());

         matrix_type MG = matrix_type::Zero(grad_basis_size, grad_basis_size);


         auto grad_quadpoints = m_bqd.grad_quadrature.integrate(msh, cl);
         for (auto& qp : grad_quadpoints)
         {
            auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());
            assert(grad_basis_size == gphi.size());

            auto dphi = m_bqd.cell_basis.eval_gradients(msh, cl, qp.point(), 0, cell_degree);
            assert(cell_basis_size == dphi.rows());
            assert(dphi.cols() == DIM);



            for(size_t j = 0; j < grad_basis_size; j++){
               for(size_t i = j; i < grad_basis_size; i++){
                  MG(i,j) += qp.weight() * mm_prod(gphi[i], gphi[j]);
               }
            }

            auto dphi_j(gphi[0]);
            for(size_t j = 0; j < cell_basis_size; j++){
               //on converti dphi_j
               dphi_j.Zero();
               for(size_t k = 0; k < DIM; k++){
                  dphi_j(k) = dphi(j,k);
               }

               for(size_t i = 0; i < grad_basis_size; i++){
                  BG(i,j) += qp.weight() * mm_prod(gphi[i], dphi_j);
               }
            }
         }// end qp

         // lower part MG
         for(size_t i = 0; i <  grad_basis_size; i++){
            for(size_t j = i; j < grad_basis_size; j++){
               MG(i,j) = MG(j,i);
            }
         }


         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto current_face_range = dsr.face_range(face_i);
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);
            auto face_quadpoints = face_quadrature.integrate(msh, fc);

            for (auto& qp : face_quadpoints)
            {
               matrix_type c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree);
               auto gphi = m_bqd.grad_basis.eval_functions(msh, cl, qp.point());

               // tau.n
               matrix_type gphi_n = matrix_type::Zero(grad_basis_size, 1);
               for(size_t i = 0; i < grad_basis_size; i++)
                     gphi_n(i,0) = mm_prod(gphi[i] ,n);

               matrix_type T = qp.weight() * gphi_n * c_phi.transpose();

               assert(T.rows() == grad_basis_size);
               assert(T.cols() == cell_basis_size);

               BG.block(0, 0, grad_basis_size, cell_basis_size) -= T;

               matrix_type f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point(), 0, face_degree);
               matrix_type F = qp.weight() * gphi_n * f_phi.transpose();

               assert(F.rows() == grad_basis_size);
               assert(F.cols() == current_face_range.size());

               BG.block(0, current_face_range.min(),
                        grad_basis_size, current_face_range.size()) += F;
            }
         }


         oper  = MG.llt().solve(BG);    // GT
         if(compute_data)
            data  = BG.transpose() * oper;  // A

            assert(oper.rows() == grad_basis_size);
         assert(oper.cols() == dsr.total_size());
      }
   };



} // namespace disk
