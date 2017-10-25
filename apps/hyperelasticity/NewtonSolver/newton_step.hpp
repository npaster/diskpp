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

 // NewtonRaphson_step

#pragma once

#include <iostream>

#include <sstream>
#include <string>
#include <fstream>


#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "hho/hho.hpp"
#include "hho/hho_vector.hpp"
#include "hho/hho_vector_bq.hpp"
#include "hho/hho_nl_vector.hpp"
#include "../ElasticityParameters.hpp"
#include "../BoundaryConditions.hpp"
#include "../Informations.hpp"
#include "../Parameters.hpp"
#include "../hyperelasticity_elementary_computation.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>



template<typename BQData>
class NewtonRaphson_step_hyperelasticity
{
   typedef typename BQData::mesh_type          mesh_type;
   typedef typename mesh_type::scalar_type     scalar_type;

   typedef dynamic_matrix<scalar_type>         matrix_dynamic;
   typedef dynamic_vector<scalar_type>         vector_dynamic;


   typedef disk::gradient_reconstruction_elas_full_bq<BQData>               gradrec_type;

   typedef disk::displacement_reconstruction_elas_bq<BQData>                deplrec_type;

   typedef Hyperelasticity::Hyperelasticity<BQData>                         hyperelasticity_type;

   typedef disk::elas_like_stabilization_PIKF_bq<BQData>                    stab_PIKF_type;
   typedef disk::elas_like_stabilization_bq<BQData>                         stab_HHO_type;
   typedef disk::diffusion_like_static_condensation_bq<BQData>              statcond_type;
   typedef disk::assembler_nl_vector_bq<BQData>                             assembler_type;


   typedef typename assembler_type::sparse_matrix_type       sparse_matrix_type;
   typedef typename assembler_type::vector_type              vector_type;

   typedef ParamRun<scalar_type>    param_type;
   #ifdef HAVE_INTEL_MKL
   typedef Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver_type;
   #else
   typedef Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver_type;
   #endif

   sparse_matrix_type     m_system_matrix;
   vector_type            m_system_rhs, m_system_solution;

   solver_type             solver;

   const BQData&                               m_bqd;

   const mesh_type& m_msh;

   std::vector<vector_dynamic>        m_solution_data, m_solution_cells, m_solution_faces, m_solution_lagr;

   std::vector<vector_dynamic>         m_RT;
   std::vector<matrix_dynamic>         m_KTT, m_KTF;
   bool m_verbose;


   ElasticityParameters m_elas_param;
   const param_type& m_rp;

   scalar_type initial_residual;

   BoundaryConditions m_boundary_condition;

public:
   NewtonRaphson_step_hyperelasticity(const mesh_type& msh, const BQData& bqd,
                                       const param_type& rp,
                                       const ElasticityParameters elas_param,
                                       const BoundaryConditions& boundary_conditions)
   : m_msh(msh), m_verbose(false), m_rp(rp), m_elas_param(elas_param), m_bqd(bqd),
     m_boundary_condition(boundary_conditions)
   {}

   bool    verbose(void) const     { return m_verbose; }
   void    verbose(bool v)         { m_verbose = v; }

   void
   initialize( const std::vector<vector_dynamic>& initial_solution_cells,
               const std::vector<vector_dynamic>& initial_solution_faces,
               const std::vector<vector_dynamic>& initial_solution_lagr,
               const std::vector<vector_dynamic>& initial_solution)
   {
      m_solution_cells.clear();
      m_solution_cells = initial_solution_cells;
      assert(m_msh.cells_size() == m_solution_cells.size());

      m_solution_faces.clear();
      m_solution_faces = initial_solution_faces;
      assert(m_msh.faces_size() == m_solution_faces.size());

      m_solution_lagr.clear();
      m_solution_lagr = initial_solution_lagr;

      m_solution_data.clear();
      m_solution_data = initial_solution;
      assert(m_msh.cells_size() == m_solution_data.size());

      m_RT.clear();
      m_RT.resize(m_msh.cells_size());

      m_KTF.clear();
      m_KTF.resize(m_msh.cells_size());

      m_KTT.clear();
      m_KTT.resize(m_msh.cells_size());
   }


//     void
//     saveMatrix(const std::string& filename)
//     {
//        std::ofstream fichier(filename, std::ios::out | std::ios::trunc);
//
//        for (int k=0; k<m_system_matrix.outerSize(); ++k)
//           for (typename Eigen::SparseMatrix<scalar_type>::InnerIterator it(m_system_matrix,k); it; ++it)
//          {
//             fichier << it.row() << " ; " << it.col() << " ; " << it.value() << std::endl;
//          }
//
//        fichier.close();
//     }
//
//     size_t
//     test_aurrichio(void)
//     {
//
//       Eigen::SelfAdjointEigenSolver<sparse_matrix_type> es;
//       es.compute(m_system_matrix);
//
//       if(es.info() != Eigen::Success) {
//          throw std::invalid_argument("ERROR: Could not compute the eigenvalues");
//       }
//
//       auto ev = es.eigenvalues();
//
//       size_t nb_negative_eigenvalue(0);
//
//       size_t i_ev = 0;
//       while(ev(i_ev) <= scalar_type(0.0))
//       {
//          nb_negative_eigenvalue++;
//          i_ev++;
//       }
//
//       const scalar_type min_ev = ev.minCoeff();
//       const scalar_type max_ev = ev.maxCoeff();
//
//       if(m_verbose){
//          std::cout << "******* Eigenvalues test ********" << std::endl;
//          std::cout << "Number of eigenvalues: " << ev.size()  << std::endl;
//          std::cout << "Number of negative eigenvalues: " << nb_negative_eigenvalue  << std::endl;
//          std::cout << "Maximum eigenvalue: " << max_ev << std::endl;
//          std::cout << "Minimum eigenvalue: " << min_ev << std::endl;
//          std::cout << "Conditionning number: " << std::abs(max_ev/min_ev) << std::endl;
//       }
//
//       return nb_negative_eigenvalue;
//     }

   SolveInfo
   solve(void)
   {
      #ifdef HAVE_INTEL_MKL
      Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
      //solver.pardisoParameterArray();[59] = 2; //out of core
      #else
      Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
      #endif

      timecounter tc;

      tc.tic();
      solver.analyzePattern(m_system_matrix);
      solver.factorize(m_system_matrix);

      if(solver.info() != Eigen::Success) {
         throw std::invalid_argument("ERROR: Could not factorize the matrix");
      }

      m_system_solution = solver.solve(m_system_rhs);
      if(solver.info() != Eigen::Success) {
         throw std::invalid_argument("ERROR: Could not solve the linear system");
      }
      tc.toc();

      return SolveInfo(m_system_matrix.rows(), m_system_matrix.nonZeros(), tc.to_double());
   }

   // Assemble
   template<typename LoadFunction, typename BoundaryConditionFunction, typename NeumannFunction>
   AssemblyInfo
   assemble(const LoadFunction& lf, const BoundaryConditionFunction& bcf, const NeumannFunction& g,
             const std::vector<matrix_dynamic>& gradient_precomputed)
   {
      // Define function
      gradrec_type gradrec(m_bqd);
      stab_HHO_type stab_HHO(m_bqd);
      hyperelasticity_type hyperelasticity(m_bqd);
      deplrec_type deplrec(m_bqd);
      stab_PIKF_type stab_PIKF(m_bqd);
      statcond_type statcond(m_bqd);
      assembler_type assembler(m_msh, m_bqd, m_boundary_condition);

      AssemblyInfo ai;

      timecounter tc, ttot;

      ttot.tic();
      size_t cell_i = 0;

      for (auto& cl : m_msh)
      {
         // Gradient Reconstruction
         matrix_dynamic GT;
         tc.tic();
         if(m_rp.m_precomputation){
            GT = gradient_precomputed[cell_i];
         }
         else{
            gradrec.compute(m_msh, cl, false);
            GT = gradrec.oper;
         }
         tc.toc();
         ai.m_time_gradrec += tc.to_double();

         // Stabilisation
         tc.tic();
         if(m_rp.m_stab){
            switch (m_rp.m_stab_type) {
               case PIKF:
               {
                  stab_PIKF.compute(m_msh, cl);

                  break;
               }
               case HHO:
               {
                  deplrec.compute(m_msh, cl);
                  stab_HHO.compute(m_msh, cl, deplrec.oper);

                  break;
               }
               case NOTHING:
               {
                  break;
               }
               default:
                  throw std::invalid_argument("Unknown stabilization");
            }
         }
         tc.toc();
         ai.m_time_stab += tc.to_double();

         // Begin Assembly
         // Build rhs and lhs

         // Mechanical Computation
         tc.tic();
         hyperelasticity.compute(m_msh, cl, lf, GT, m_solution_data.at(cell_i), m_elas_param);
         dynamic_matrix<scalar_type> lhs = hyperelasticity.K_int;
         dynamic_vector<scalar_type> rhs = hyperelasticity.RTF;

         tc.toc();
         ai.m_time_elem += tc.to_double();
         ai.m_time_law += hyperelasticity.time_law;

         // Stabilisation Contribution
         tc.tic();
         if(m_rp.m_stab){
            switch (m_rp.m_stab_type) {
               case PIKF:
               {
                  assert( hyperelasticity.K_int.rows() == stab_PIKF.data.rows());
                  assert( hyperelasticity.K_int.cols() == stab_PIKF.data.cols());
                  assert( hyperelasticity.RTF.rows() == (stab_PIKF.data * m_solution_data.at(cell_i)).rows());
                  assert( hyperelasticity.RTF.cols() == (stab_PIKF.data * m_solution_data.at(cell_i)).cols());

                  lhs += m_rp.m_beta * stab_PIKF.data;
                  rhs -= m_rp.m_beta * stab_PIKF.data * m_solution_data.at(cell_i);
                  break;
               }
               case HHO:
               {
                  assert( hyperelasticity.K_int.rows() == stab_HHO.data.rows());
                  assert( hyperelasticity.K_int.cols() == stab_HHO.data.cols());
                  assert( hyperelasticity.RTF.rows() == (stab_HHO.data * m_solution_data.at(cell_i)).rows());
                  assert( hyperelasticity.RTF.cols() == (stab_HHO.data * m_solution_data.at(cell_i)).cols());

                  lhs += m_rp.m_beta * stab_HHO.data;
                  rhs -= m_rp.m_beta * stab_HHO.data * m_solution_data.at(cell_i);
                  break;
               }
               case NOTHING:
               {
                  break;
               }
               default:
                  throw std::invalid_argument("Unknown stabilization");
            }
         }
         tc.toc();
         ai.m_time_stab += tc.to_double();

         // Static Condensation
         tc.tic();
         auto scnp = statcond.compute(m_msh, cl, lhs, rhs, true);
         m_KTT[cell_i] = statcond.KTT;
         m_KTF[cell_i] = statcond.KTF;
         m_RT[cell_i] = statcond.RT;
         tc.toc();
         ai.m_time_statcond += tc.to_double();

         assembler.assemble(m_msh, cl, scnp);

         cell_i++;
      }

      // Impose Boundary Conditions
      assembler.impose_boundary_conditions(m_msh, bcf, g, m_solution_faces, m_solution_lagr, m_boundary_condition);

      // Assemble the global system
      assembler.finalize(m_system_matrix, m_system_rhs);

      // Update Infos
      ttot.toc();
      ai.m_time_assembly = ttot.to_double();
      ai.m_linear_system_size = m_system_matrix.rows();
      return ai;
   }

   // Post-processing
   void postprocess()
   {
      // Number of Unknowns by cell and face
      const size_t cbs = (m_bqd.cell_basis.range(0, m_bqd.cell_degree())).size();
      const size_t fbs = (m_bqd.face_basis.range(0, m_bqd.face_degree())).size();

      // Update  unknowns
      // Update face Uf^{i+1} = Uf^i + delta Uf^i
      for(size_t i = 0; i < m_solution_faces.size(); i++){
         assert(m_solution_faces.at(i).size() == fbs);
         m_solution_faces.at(i) += m_system_solution.block(i * fbs, 0, fbs, 1);
      }

      // Update lagrangian L^{i+1} = L^i + delta L^i
      const size_t lagrange_offset = m_solution_faces.size() * fbs;
      const size_t num_lagr_dofs = fbs/m_msh.dimension;

      for(size_t i = 0; i < m_solution_lagr.size(); i++){
         const size_t pos = lagrange_offset + m_boundary_condition.begin_lag_conditions_faceI(i);
         const size_t size = num_lagr_dofs * m_boundary_condition.nb_lag_conditions_faceI(i);
         m_solution_lagr.at(i) += m_system_solution.block(lagrange_offset + i * size, 0, size, 1);
      }

      // Update cell
      size_t cell_i(0);
      for (auto& cl : m_msh)
      {
         // Extract the solution
         const auto fcs = faces(m_msh, cl);
         const auto num_faces = fcs.size();

         dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(num_faces*fbs);

         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            const auto fc = fcs[face_i];
            const auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), fc);
            if (!eid.first)
               throw std::invalid_argument("This is a bug: face not found");

            const auto face_id = eid.second;

            xFs.block(face_i * fbs, 0, fbs, 1) = m_system_solution.block(face_id * fbs, 0, fbs, 1);
         }

         auto K_TT_ldlt = m_KTT[cell_i].llt();
         const dynamic_vector<scalar_type> xT = - K_TT_ldlt.solve(-m_RT[cell_i] + m_KTF[cell_i] * xFs);

         assert(xT.size() == cbs);
         assert(m_solution_data.at(cell_i).size() == xT.size() + xFs.size());
         // Update element U^{i+1} = U^i + delta U^i ///
         (m_solution_data.at(cell_i)).block(0,0,cbs,1) += xT;
         (m_solution_data.at(cell_i)).block(cbs,0,fbs,1) += xFs;

         // Update Cell Uc^{i+1} = Uc^i + delta Uc^i ///
         assert(m_solution_cells.at(cell_i).size() == cbs);
         m_solution_cells.at(cell_i) += xT;

         cell_i++;
      }
   }


   std::array<scalar_type, 2>
   compute_energy() const
   {
      gradrec_type gradrec(m_bqd);
      deplrec_type deplrec(m_bqd);
      stab_HHO_type stab_HHO(m_bqd);
      stab_PIKF_type stab_PIKF(m_bqd);

      std::array<scalar_type, 2> ret = {0.0, 0.0};
      scalar_type energy_stab(0.0);

      const NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);

      size_t i = 0;
      for (auto& cl : m_msh)
      {
         gradrec.compute(m_msh, cl, false);

         // Energie in the stabilisation
         if(m_rp.m_stab){
            switch (m_rp.m_stab_type) {
               case PIKF:
               {
                  stab_PIKF.compute(m_msh, cl);
                  energy_stab = m_rp.m_beta *  m_solution_data.at(i).dot(stab_PIKF.data * m_solution_data.at(i));
                  break;
               }
               case HHO:
               {
                  deplrec.compute(m_msh, cl);
                  stab_HHO.compute(m_msh, cl, deplrec.oper);
                  energy_stab = m_rp.m_beta *  m_solution_data.at(i).dot(stab_HHO.data * m_solution_data.at(i));
                  break;
               }
               case NOTHING:
               {
                  break;
               }
               default:
                  std::cout << "Unknown Stabilisation " << m_rp.m_stab_type << std::endl;
                  throw std::invalid_argument("Unknown stabilization");
            }
         }

         const vector_type GT_uTF = gradrec.oper * m_solution_data.at(i);
         const auto grad_quadpoints = m_bqd.grad_quadrature.integrate(m_msh, cl);

         for (auto& qp : grad_quadpoints)
         {
            const auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, qp.point());

            // Compute local gradient and norm
            const auto GT_iqn = disk::compute_gradient_matrix_pt(GT_uTF, gphi);
            const auto FT_iqn = compute_FTensor(GT_iqn);

            const scalar_type energy = law.compute_energy(FT_iqn);

            //Elastic energy
            ret[0] += qp.weight() * energy;

            // Energie in the stabilisation
            if(m_rp.m_stab){
               ret[1] += qp.weight() * energy_stab;
            }
         }
         i++;
      }

      return ret;
   }


   bool test_convergence(const scalar_type epsilon, const size_t iter, scalar_type& error)
   {
      if(iter == 0){
         initial_residual = m_system_rhs.norm();
      }

      const scalar_type residual = m_system_rhs.norm();
      scalar_type max_error(0.0);
      for(size_t i = 0; i < m_system_rhs.size(); i++)
         max_error = std::max( max_error, std::abs(m_system_rhs(i)));

      scalar_type relative_error(1.E4);

      if(initial_residual == static_cast<scalar_type>(0))
      {
         relative_error = static_cast<scalar_type>(0);
         max_error = static_cast<scalar_type>(0);
      }
      else {
         relative_error = residual / initial_residual;
      }

      const scalar_type error_incr = m_system_solution.norm();

      scalar_type error_un(0.0);

      for (size_t i = 0; i < m_solution_faces.size(); i++) {
         scalar_type norm = m_solution_faces[i].norm();
         error_un += norm * norm;
      }

      for (size_t i = 0; i < m_solution_lagr.size(); i++) {
         scalar_type norm = m_solution_lagr[i].norm();
         error_un += norm * norm;
      }

      error_un = std::sqrt(error_un);

      if(error_un == static_cast<scalar_type>(0)){
         error_un = static_cast<scalar_type>(10E6);
      }

      scalar_type relative_displ = error_incr/error_un;

      if(iter == 0){
         relative_displ = static_cast<scalar_type>(1.0);
      }

      const size_t nb_faces_dof = m_bqd.face_basis.size() * m_msh.faces_size();

      const scalar_type norm_depl = (m_system_rhs.head(nb_faces_dof)).norm();
      const scalar_type norm_lag = (m_system_rhs.tail(m_system_rhs.size() - nb_faces_dof)).norm();

      if(m_verbose){
         std::string s_iter = "   " + std::to_string(iter) + "               ";
         s_iter.resize(9);

         if(iter == 0){
            std::cout << "--------------------------------------------------------------------------------------------------------------------------------" << std::endl;
            std::cout << "| Iteration | Norme l2 incr | Relative depl |  Residual l2  | Relative error | Maximum error | Residual face  |  Residual BC   |" << std::endl;
            std::cout << "--------------------------------------------------------------------------------------------------------------------------------" << std::endl;

         }
         std::ios::fmtflags f( std::cout.flags() );
         std::cout.precision(5);
         std::cout.setf(std::iostream::scientific, std::iostream::floatfield);
         std::cout << "| " << s_iter << " |   " << error_incr << " |   " << relative_displ << " |   " << residual << " |   " << relative_error << "  |  " << max_error << "  |   "
         << norm_depl << "  |   " << norm_lag << "  |" << std::endl;
         std::cout << "--------------------------------------------------------------------------------------------------------------------------------" << std::endl;
         std::cout.flags( f );
      }

      error = std::min(relative_displ, relative_error);

      if(error <= epsilon ){
         return true;
      }
      else {
         return false;
      }

   }

   void
   save_solutions( std::vector<vector_dynamic>& solution_cells,
                   std::vector<vector_dynamic>& solution_faces,
                   std::vector<vector_dynamic>& solution_lagr,
                   std::vector<vector_dynamic>& solution) const
    {
      solution_cells.clear();
      solution_cells = m_solution_cells;
      assert(m_solution_cells.size() == solution_cells.size());

      solution_faces.clear();
      solution_faces = m_solution_faces;
      assert(m_solution_faces.size() == solution_faces.size());

      solution_lagr.clear();
      solution_lagr = m_solution_lagr;
      assert(m_solution_lagr.size() == solution_lagr.size());

      solution.clear();
      solution = m_solution_data;
      assert(m_solution_data.size() == solution.size());
    }
};
