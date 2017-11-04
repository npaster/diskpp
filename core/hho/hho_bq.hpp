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

namespace hho{

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
      face_quad_type      face_max_quadrature;

   private:
      size_t  m_cell_degree, m_face_degree;

      void init(void)
      {
         cell_basis          = cell_basis_type(m_cell_degree+1);
         face_basis          = face_basis_type(m_face_degree);
         cell_quadrature     = cell_quad_type(2*(m_cell_degree+1));
         face_quadrature     = face_quad_type(2*m_face_degree);
         face_max_quadrature = face_quad_type(2* std::max(m_face_degree,m_cell_degree+1));
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

      cell_basis_type     cell_basis;
      face_basis_type     face_basis;
      cell_quad_type      cell_quadrature;
      face_quad_type      face_quadrature;
      grad_basis_type     grad_basis;
      cell_quad_type      grad_quadrature;

      cell_quad_type      grad_cell_quadrature;
      face_quad_type      grad_face_quadrature;
      face_quad_type      grad_face_max_quadrature;

   private:
      size_t  m_cell_degree, m_face_degree, m_grad_degree;

      void init(void)
      {
         cell_basis          = cell_basis_type(m_cell_degree);
         face_basis          = face_basis_type(m_face_degree);
         grad_basis          = grad_basis_type(m_grad_degree);
         cell_quadrature     = cell_quad_type(2 * m_cell_degree);
         face_quadrature     = face_quad_type(2 * m_face_degree);
         grad_quadrature     = cell_quad_type(2 * m_grad_degree);
         grad_cell_quadrature     = cell_quad_type(m_grad_degree + m_cell_degree);
         grad_face_quadrature     = face_quad_type(m_grad_degree + m_face_degree);
         grad_face_max_quadrature     = face_quad_type(m_grad_degree + std::max(m_face_degree, m_cell_degree));
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



   template<typename BQData, typename Function>
   dynamic_vector<typename BQData::mesh_type::scalar_type>
   compute_rhs_bq(const typename BQData::mesh_type& msh,
                  const typename BQData::mesh_type::cell& cl,
                  const Function& f, const BQData& bqd)
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef dynamic_vector<scalar_type> vector_type;

      vector_type ret = vector_type::Zero(bqd.cell_basis.size());

      const auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
      for (auto& qp : cell_quadpoints)
      {
         auto phi = bqd.cell_basis.eval_functions(msh, cl, qp.point());
         auto fval = f(qp.point());
         for (size_t i = 0; i < bqd.cell_basis.size(); i++)
            ret(i) += qp.weight() * mm_prod(fval, phi[i]);
      }

      return ret;
   }


   template<typename BQData>
   class assembler_bq
   {
      typedef typename BQData::mesh_type          mesh_type;
      typedef typename mesh_type::scalar_type     scalar_type;
      typedef typename mesh_type::cell            cell_type;
      typedef typename mesh_type::face            face_type;

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

      assembler_bq()                 = delete;

      assembler_bq(const mesh_type& msh, const BQData& bqd)
      : m_bqd(bqd)
      {
         m_num_unknowns = m_bqd.face_basis.size() * (msh.faces_size() + msh.boundary_faces_size());
         matrix = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
         rhs = vector_type::Zero(m_num_unknowns);
      }

      template<typename LocalContrib>
      void
      assemble(const mesh_type& msh, const cell_type& cl, const LocalContrib& lc)
      {
         const size_t face_basis_size = m_bqd.face_basis.size();

         const auto fcs = faces(msh, cl);
         std::vector<size_t> l2g(fcs.size() * face_basis_size);

         assert(face_basis_size == face_basis_size);

         for (size_t face_i = 0; face_i < fcs.size(); face_i++)
         {
            const auto fc = fcs[face_i];
            const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
               throw std::invalid_argument("This is a bug: face not found");

            const auto face_id = eid.second;

            const auto face_offset = face_id * face_basis_size;

            const auto pos = face_i * face_basis_size;

            for (size_t i = 0; i < face_basis_size; i++)
               l2g[pos+i] = face_offset+i;
         }

         assert(lc.first.rows() == fcs.size()*face_basis_size);
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
         const size_t fbs = m_bqd.face_basis.size();
         size_t face_i = 0;
         for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
         {
            auto bfc = *itor;

            const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first)
               throw std::invalid_argument("This is a bug: face not found");

            const auto face_id = eid.second;

            const auto face_offset = face_id * fbs;
            const auto face_offset_lagrange = (msh.faces_size() + face_i) * fbs;

            const auto fqd = m_bqd.face_quadrature.integrate(msh, bfc);

            matrix_type MFF     = matrix_type::Zero(fbs, fbs);
            vector_type rhs_f   = vector_type::Zero(fbs);

            for (auto& qp : fqd)
            {
               const auto f_phi = m_bqd.face_basis.eval_functions(msh, bfc, qp.point());
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

      template<typename Function>
      void
      impose_boundary_conditions_nl(const mesh_type& msh, const Function& bc, const std::vector<vector_type>& sol_faces,
                                 const std::vector<vector_type>& sol_lagr)
      {
         const size_t fbs = m_bqd.face_basis.size();
         size_t face_i = 0;
         for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
         {
            auto bfc = *itor;

            const auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first)
               throw std::invalid_argument("This is a bug: face not found");

            const auto face_id = eid.second;

            const auto face_offset = face_id * fbs;
            const auto face_offset_lagrange = (msh.faces_size() + face_i) * fbs;

            const auto fqd = m_bqd.face_quadrature.integrate(msh, bfc);

            matrix_type MFF     = matrix_type::Zero(fbs, fbs);
            vector_type rhs_f   = vector_type::Zero(fbs);

            for (auto& qp : fqd)
            {
               const auto f_phi = m_bqd.face_basis.eval_functions(msh, bfc, qp.point());
               MFF += qp.weight() * f_phi * f_phi.transpose();
               rhs_f += qp.weight() * f_phi * bc(qp.point());
            }


            rhs_f -= MFF * sol_faces.at(face_id);

            const vector_type rhs_l = MFF * sol_lagr.at(face_i);

            for (size_t i = 0; i < fbs; i++)
               rhs(face_offset+i) -= rhs_l(i);

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

} // namespace hho

} // namespace disk
