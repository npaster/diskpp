/*
*       /\        Matteo Cicuttin (C) 2016, 2017
*      /__\       matteo.cicuttin@enpc.fr
*     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
*    /\    /\
*   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
*  /_\/_\/_\/_\   methods.
*
* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, You can obtain one at http://mozilla.org/MPL/2.0/.
*
* If you use this code or parts of it for scientific publications, you
* are required to cite it as following:
*
* Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
* polytopal meshes using generic programming.
* M. Cicuttin, D. A. Di Pietro, A. Ern.
* Journal of Computational and Applied Mathematics.
* DOI: 10.1016/j.cam.2017.09.017
*/

#pragma once

#include "common/eigen.hpp"
#include "hho/hho_bq.hpp"

namespace disk{

   namespace hho{

      template<typename BQData>
      struct hho_stabilization_bq
      {
          static_assert(sizeof(BQData) == -1, "hho_stabilization: not suitable for the requested kind of BQData");
      };

      //partial specialization for disk::scaled_monomial_scalar_basis
      template<typename MeshType, template<typename, typename> class QuadratureType>
      class hho_stabilization_bq<basis_quadrature_data<  MeshType,
                                                         disk::scaled_monomial_scalar_basis,
                                                         QuadratureType> >
      {
      private:
         typedef basis_quadrature_data<   MeshType,
                                          disk::scaled_monomial_scalar_basis,
                                          QuadratureType>      BQData_type;

         typedef typename BQData_type::mesh_type               mesh_type;
         typedef typename mesh_type::scalar_type               scalar_type;
         typedef typename mesh_type::cell                      cell_type;

         typedef dynamic_matrix<scalar_type>                   matrix_type;
         typedef dynamic_vector<scalar_type>                   vector_type;

         const BQData_type&                                    m_bqd;

      public:
         matrix_type     data;

         hho_stabilization_bq(const BQData_type& bqd) : m_bqd(bqd)
         {}

         void compute(const mesh_type& msh, const cell_type& cl, const matrix_type& gradrec_oper)
         {
            const auto cell_degree = m_bqd.cell_degree();
            const auto face_degree = m_bqd.face_degree();
            const auto cell_basis_size = m_bqd.cell_basis.computed_size();
            const auto face_basis_size = howmany_dofs(m_bqd.face_basis);

            matrix_type mass_mat = matrix_type::Zero(cell_basis_size, cell_basis_size);

            const auto cell_quadpoints = m_bqd.cell_quadrature.integrate(msh, cl);
            for (auto& qp : cell_quadpoints)
            {
               const auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree+1);
               assert(c_phi.rows() == cell_basis_size);
               mass_mat += qp.weight() * c_phi * c_phi.transpose();
            }

            const auto zero_range         = m_bqd.cell_basis.range(0, cell_degree);
            const auto one_range          = m_bqd.cell_basis.range(1, cell_degree+1);

            // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

            //Step 1: compute \pi_T^k p_T^k v (third term).
            const matrix_type M1 = take(mass_mat, zero_range, zero_range);
            const matrix_type M2 = take(mass_mat, zero_range, one_range);
            matrix_type proj1 = -M1.llt().solve(M2*gradrec_oper);

            //Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
            const matrix_type I_T = matrix_type::Identity(zero_range.size(), zero_range.size());
            proj1.block(0, 0, zero_range.size(), zero_range.size()) += I_T;

            const auto fcs = faces(msh, cl);
            const auto num_faces = fcs.size();
            const auto num_cell_dofs = howmany_dofs(m_bqd.cell_basis);
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

               const auto face_quadpoints = m_bqd.face_max_quadrature.integrate(msh, fc);
               for (auto& qp : face_quadpoints)
               {
                  const auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());
                  const auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point(), 0, cell_degree+1);
                  assert(c_phi.size() == cell_basis_size);
                  assert(f_phi.size() == face_basis_size);
                  const auto q_f_phi = qp.weight() * f_phi;

                  face_mass_matrix += q_f_phi * f_phi.transpose();
                  face_trace_matrix += q_f_phi * c_phi.transpose();
               }

               Eigen::LLT<matrix_type> piKF;
               piKF.compute(face_mass_matrix);

               // Step 3a: \pi_F^k( v_F - p_T^k v )
               const auto face_range = current_face_range.remove_offset();
               const matrix_type MR1 = take(face_trace_matrix, face_range, one_range);
               matrix_type proj2 = piKF.solve(MR1*gradrec_oper);
               const matrix_type I_F = matrix_type::Identity(face_basis_size, face_basis_size);
               const auto block_offset = current_face_range.min();
               proj2.block(0, block_offset, face_basis_size, face_basis_size) -= I_F;

               // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
               const matrix_type MR2 = take(face_trace_matrix, face_range, zero_range);
               const matrix_type proj3 = piKF.solve(MR2*proj1);

               const matrix_type BRF = proj2 + proj3;

               data += BRF.transpose() * face_mass_matrix * BRF / h;
            }
         }
      };


      //partial specialization for disk::scaled_monomial_scalar_basis
      template<typename MeshType, template<typename, typename> class GradBasisType,
      template<typename, typename> class QuadratureType>
      class hho_stabilization_bq<basis_quadrature_data_full<MeshType,
                                                            disk::scaled_monomial_scalar_basis,
                                                            GradBasisType,
                                                            QuadratureType> >
      {
      private:
         typedef basis_quadrature_data_full< MeshType,
                                             disk::scaled_monomial_scalar_basis,
                                             GradBasisType,
                                             QuadratureType>   BQData_type;

         typedef basis_quadrature_data<MeshType,
                                       disk::scaled_monomial_scalar_basis,
                                       QuadratureType>         BQData_type2;

         typedef typename BQData_type::mesh_type               mesh_type;
         typedef typename mesh_type::scalar_type               scalar_type;
         typedef typename mesh_type::cell                      cell_type;

         typedef dynamic_matrix<scalar_type>                   matrix_type;
         typedef dynamic_vector<scalar_type>                   vector_type;

         const BQData_type&                                    m_bqd;
         BQData_type2                                          m_bqd2;

      public:
         matrix_type     data;

         hho_stabilization_bq(const BQData_type& bqd) : m_bqd(bqd)
         {
            m_bqd2 = BQData_type2(bqd.cell_degree()+1, bqd.face_degree());
         }

         void compute(const mesh_type& msh, const cell_type& cl, const matrix_type& gradrec_oper)
         {
            typedef hho_stabilization_bq<BQData_type2> hho_stabtype;

            hho_stabtype stabilization(m_bqd2);
            stabilization.compute(msh, cl, gradrec_oper);
            data = stabilization.data;
         }
      };


      //partial specialization for disk::scaled_monomial_vector_basis
      template<typename MeshType, template<typename, typename> class QuadratureType>
      class hho_stabilization_bq<basis_quadrature_data<  MeshType,
                                                         disk::scaled_monomial_vector_basis,
                                                         QuadratureType> >
      {
      private:
         typedef basis_quadrature_data<   MeshType,
                                          disk::scaled_monomial_vector_basis,
                                          QuadratureType>      BQData_type;

         typedef typename BQData_type::mesh_type               mesh_type;
         typedef typename mesh_type::scalar_type               scalar_type;
         typedef typename mesh_type::cell                      cell_type;

         typedef dynamic_matrix<scalar_type>                   matrix_type;
         typedef dynamic_vector<scalar_type>                   vector_type;

         const BQData_type&                                    m_bqd;

      public:
         matrix_type     data;

         hho_stabilization_bq(const BQData_type& bqd) : m_bqd(bqd)
         {}

         void compute(const mesh_type& msh, const cell_type& cl, const matrix_type& gradrec_oper)
         {
            const size_t DIM = msh.dimension;
            const auto cell_degree = m_bqd.cell_degree();
            const auto face_degree = m_bqd.face_degree();
            const auto cell_basis_size = m_bqd.cell_basis.range(0,cell_degree+1).size();
            const auto face_basis_size = m_bqd.face_basis.range(0,face_degree).size();

            matrix_type mass_mat = matrix_type::Zero(cell_basis_size, cell_basis_size);

            const auto cell_quadpoints = m_bqd.cell_quadrature.integrate(msh, cl);
            for (auto& qp : cell_quadpoints)
            {
               const auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point());
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
std::cout << "mass fin" << '\n';
            const auto zero_range         = m_bqd.cell_basis.range(0, cell_degree);
            const auto one_range          = m_bqd.cell_basis.range(1, cell_degree+1);

            // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

            //Step 1: compute \pi_T^k p_T^k v (third term).
            const matrix_type M1 = take(mass_mat, zero_range, zero_range);
            const matrix_type M2 = take(mass_mat, zero_range, one_range);
            std::cout << "3" << '\n';
            std::cout << M2.rows() <<" " << M2.cols() << '\n';
            std::cout << gradrec_oper.rows() <<" " << gradrec_oper.cols() << '\n';
            matrix_type proj1 = -M1.llt().solve(M2*gradrec_oper);
            std::cout << "2" << '\n';

            //Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
            const matrix_type I_T = matrix_type::Identity(zero_range.size(), zero_range.size());
            proj1.block(0, 0, zero_range.size(), zero_range.size()) += I_T;

            const auto fcs = faces(msh, cl);
            const auto num_faces = fcs.size();
            const size_t num_cell_dofs = m_bqd.cell_basis.range(0, cell_degree).size();

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
std::cout << "face" << '\n';
               const auto face_quadpoints = m_bqd.face_max_quadrature.integrate(msh, fc);
               for (auto& qp : face_quadpoints)
               {
                  const auto f_phi = m_bqd.face_basis.eval_functions(msh, fc, qp.point());
                  const auto c_phi = m_bqd.cell_basis.eval_functions(msh, cl, qp.point());

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
std::cout << "é" << '\n';
                     for(size_t j = 0; j < cell_basis_size; j++){
                        for(size_t i = 0; i< face_basis_size; i++){
                           face_trace_matrix(i,j) += qp.weight() * mm_prod(f_phi[i], c_phi[j]);
                        }
                     }
               }

               //lower part
               for (size_t i = 0; i < face_basis_size; i++)
                  for (size_t j = i; j < face_basis_size; j++)
                     face_mass_matrix(i,j) = face_mass_matrix(j,i);
std::cout << "face fin" << '\n';
               Eigen::LLT<matrix_type> piKF;
               piKF.compute(face_mass_matrix);

               // Step 3a: \pi_F^k( v_F - p_T^k v )
               const auto face_range = current_face_range.remove_offset();
               const matrix_type MR1 = take(face_trace_matrix, face_range, one_range);
               matrix_type proj2 = piKF.solve(MR1*gradrec_oper);
               const matrix_type I_F = matrix_type::Identity(face_basis_size, face_basis_size);
               const auto block_offset = current_face_range.min();
               proj2.block(0, block_offset, face_basis_size, face_basis_size) -= I_F;

               // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
               const matrix_type MR2 = take(face_trace_matrix, face_range, zero_range);
               const matrix_type proj3 = piKF.solve(MR2*proj1);

               const matrix_type BRF = proj2 + proj3;

               data += BRF.transpose() * face_mass_matrix * BRF / h;
            }
         }
      };


      //partial specialization for disk::scaled_monomial_vector_basis
      template<typename MeshType, template<typename, typename> class GradBasisType,
      template<typename, typename> class QuadratureType>
      class hho_stabilization_bq<basis_quadrature_data_full<MeshType,
                                                            disk::scaled_monomial_vector_basis,
                                                            GradBasisType,
                                                            QuadratureType> >
      {
      private:
         typedef basis_quadrature_data_full< MeshType,
                                             disk::scaled_monomial_vector_basis,
                                             GradBasisType,
                                             QuadratureType>   BQData_type;

         typedef basis_quadrature_data<MeshType,
                                       disk::scaled_monomial_vector_basis,
                                       QuadratureType>         BQData_type2;

         typedef typename BQData_type::mesh_type              mesh_type;
         typedef typename mesh_type::scalar_type               scalar_type;
         typedef typename mesh_type::cell                      cell_type;

         typedef dynamic_matrix<scalar_type>                   matrix_type;
         typedef dynamic_vector<scalar_type>                   vector_type;

         const BQData_type&                                    m_bqd;
         BQData_type2                                          m_bqd2;

      public:
         matrix_type     data;

         hho_stabilization_bq(const BQData_type& bqd) : m_bqd(bqd)
         {
            m_bqd2 = BQData_type2(bqd.cell_degree()+1, bqd.face_degree());
         }

         void compute(const mesh_type& msh, const cell_type& cl, const matrix_type& gradrec_oper)
         {
            typedef hho_stabilization_bq<BQData_type2> hho_stabtype;

            hho_stabtype stabilization(m_bqd2);
            stabilization.compute(msh, cl, gradrec_oper);
            data = stabilization.data;
         }
      };

   }// hho
}//disk
