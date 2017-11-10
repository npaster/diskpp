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

#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>

#include "config.h"

#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "NewtonSolver/newton_solver.hpp"
#include "hyperelasticity_elementary_computation.hpp"
#include "ElasticityParameters.hpp"
#include "Informations.hpp"
#include "Parameters.hpp"

#include "hho/hho.hpp"
#include "hho/hho_bq.hpp"
#include "hho/hho_vector_bq.hpp"
#include "hho/hho_utils.hpp"
#include "bases/bases_utils.hpp"

#include "mechanics/stress_tensors.hpp"
#include "mechanics/deformation_tensors.hpp"
#include "mechanics/BoundaryConditions.hpp"

#include "output/gmshDisk.hpp"
#include "output/gmshConvertMesh.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>



struct time_step
{
   double time;
   size_t level;
};

template<template<typename, size_t , typename> class Mesh,
typename T, size_t DIM, typename Storage, typename Point>
class hyperelasticity_solver
{
   typedef Mesh<T, DIM, Storage>                      mesh_type;
   typedef typename mesh_type::scalar_type            scalar_type;
   typedef ParamRun<T>                                param_type;

   typedef disk::hho::basis_quadrature_data_full<mesh_type, disk::scaled_monomial_vector_basis,
   //disk::Raviart_Thomas_matrix_basis,
   disk::scaled_monomial_matrix_basis,
                                                      disk::quadrature> bqdata_type;

   typedef dynamic_matrix<scalar_type>         matrix_dynamic;
   typedef dynamic_vector<scalar_type>         vector_dynamic;

   typedef disk::hho::gradient_reconstruction_vector_full_bq<bqdata_type>     gradrec_type;
   typedef disk::hho::projector_vector_bq<bqdata_type>                        projector_type;
   typedef disk::hho::displacement_reconstruction_elas_bq<bqdata_type>           deplrec_type;

   bqdata_type     m_bqd;

   std::vector<matrix_dynamic>  m_gradient_precomputed;

   const mesh_type& m_msh;

   std::vector<vector_dynamic>         m_solution_data;
   std::vector<vector_dynamic>         m_solution_cells, m_solution_faces, m_solution_lagr;

   bool m_verbose, m_convergence;

   ElasticityParameters m_elas_param;
   param_type m_rp;

   BoundaryConditions m_boundary_condition;

   size_t total_dof_depl_static;


   void
   pre_computation()
   {
      gradrec_type gradrec(m_bqd);

      m_gradient_precomputed.clear();
      m_gradient_precomputed.reserve(m_msh.cells_size());

      for (auto& cl : m_msh)
      {
         /////// Gradient Reconstruction /////////
         gradrec.compute(m_msh, cl, false);
         m_gradient_precomputed.push_back(gradrec.oper);
      }

   }

public:
   hyperelasticity_solver(const mesh_type& msh, const param_type& rp, const ElasticityParameters elas_param,
                         const std::vector<BoundaryType>& boundary_neumann,
                         const std::vector<BoundaryType>& boundary_dirichlet)
   : m_msh(msh), m_verbose(rp.m_verbose), m_convergence(false), m_elas_param(elas_param), m_rp(rp)
   {
      int face_degree = rp.m_face_degree;
      if(face_degree < 0)
      {
         std::cout << "'face_degree' should be > 0. Reverting to 1." << std::endl;
         face_degree = 1;
      }

      m_rp.m_face_degree = face_degree;

      int cell_degree = rp.m_cell_degree;
      if( ((face_degree-1)  > cell_degree) or (cell_degree > (face_degree +1)) or (cell_degree < 0))
      {
         std::cout << "face degree: " << face_degree << std::endl;
         std::cout << "cell degree: " << cell_degree << std::endl;
         std::cout << "'cell_degree' should be 'face_degree + 1' => 'cell_degree' => 'face_degree -1'. Reverting to 'face_degree'." << std::endl;
         cell_degree = face_degree;
      }

      m_rp.m_cell_degree = cell_degree;

      int grad_degree = rp.m_grad_degree;
      if(grad_degree  < cell_degree)
      {
         std::cout << "'grad_degree' should be > 'cell_degree'. Reverting to 'cell_degree'." << std::endl;
         grad_degree = cell_degree;
      }

      m_rp.m_grad_degree = grad_degree;

      m_bqd = bqdata_type(m_rp.m_face_degree, m_rp.m_cell_degree, m_rp.m_grad_degree);

      if(m_verbose){
         m_bqd.info_degree();
         m_rp.infos();
      }

      m_boundary_condition = BoundaryConditions(msh, boundary_neumann, boundary_dirichlet);

   }

   bool    verbose(void) const     { return m_verbose; }
   void    verbose(bool v)         { m_verbose = v; }


   void
   compute_initial_state()
   {
      m_solution_data.clear();
      m_solution_cells.clear();
      m_solution_faces.clear();
      m_solution_lagr.clear();

      const size_t nb_faces_dirichlet = m_boundary_condition.nb_faces_dirichlet();
      const size_t nb_lag_conditions = m_boundary_condition.nb_lags();

      m_solution_data.reserve(m_msh.cells_size());
      m_solution_cells.reserve(m_msh.cells_size());
      m_solution_faces.reserve(m_msh.faces_size());
      m_solution_lagr.reserve(nb_faces_dirichlet);


      const size_t num_cell_dofs = (m_bqd.cell_basis.range(0, m_bqd.cell_degree())).size();
      const size_t num_face_dofs = m_bqd.face_basis.size();
      const size_t num_lagr_dofs = num_face_dofs/DIM;
      const size_t total_dof = m_msh.cells_size() * num_cell_dofs + m_msh.faces_size() * num_face_dofs;
      const size_t total_lagr = nb_lag_conditions * num_face_dofs/m_msh.dimension;

      for(auto& cl : m_msh)
      {
         const auto fcs = faces(m_msh, cl);
         const size_t num_faces = fcs.size();
         m_solution_data.push_back(vector_dynamic::Zero(num_cell_dofs + num_faces * num_face_dofs));
         m_solution_cells.push_back(vector_dynamic::Zero(num_cell_dofs));
      }

      for(size_t i = 0; i < m_msh.faces_size(); i++)
      {
         m_solution_faces.push_back(vector_dynamic::Zero(num_face_dofs));
      }

      for(size_t i = 0; i < nb_faces_dirichlet; i++){
         // std::cout << "nb lag" << m_boundary_condition.nb_lag_conditions_faceI(i)  << '\n';
         // std::cout << "taille " << num_face_dofs << " vs " << num_lagr_dofs * m_boundary_condition.nb_lag_conditions_faceI(i)  << '\n';
         m_solution_lagr.push_back(vector_dynamic::Zero(num_lagr_dofs * m_boundary_condition.nb_lag_conditions_faceI(i)));
      }

      total_dof_depl_static = m_msh.faces_size() * num_face_dofs;

      if(m_verbose){
         std::cout << "** Numbers of cells: " << m_msh.cells_size()  << std::endl;
         std::cout << "** Numbers of faces: " << m_msh.faces_size() << " ( boundary faces: "
         <<  m_msh.boundary_faces_size() << " )"  << std::endl;
         std::cout << "** Numbers of dofs: " << total_dof + total_lagr  << std::endl;
         std::cout << "** including " << total_lagr << " Lagrange multipliers"  << std::endl;
         std::cout << "** After static condensation: "  << std::endl;
         std::cout << "** Numbers of dofs: " << total_dof_depl_static + total_lagr  << std::endl;
         std::cout << "** including " << total_lagr << " Lagrange multipliers"  << std::endl;
      }

//       //provisoire
//
//       for(size_t i = 0; i < m_msh.cells_size(); i++)
//       {
//          m_solution_data[i].setConstant(1E-8);
//          m_solution_cells[i].setConstant(0.01);
//       }
// //
//       for(size_t i = 0; i < m_msh.faces_size(); i++)
//       {
//          m_solution_faces[i].setConstant(1E-8);
//       }
//
//       for(size_t i = 0; i < m_msh.boundary_faces_size(); i++){
//          m_solution_lagr[i].setConstant(1.0);
//      }
   }

   template<typename LoadFunction, typename BoundaryConditionFunction, typename NeumannFunction>
   SolverInfo
   compute(const LoadFunction& lf, const BoundaryConditionFunction& bcf, const NeumannFunction& g)
   {
      SolverInfo si;
      timecounter ttot;
      ttot.tic();

      // Initialize Newton solver
      NewtonRaphson_solver_hyperelasticity<bqdata_type>
         newton_solver(m_msh, m_bqd, m_rp, m_elas_param, m_boundary_condition);

      newton_solver.initialize(m_solution_cells, m_solution_faces,
                                 m_solution_lagr, m_solution_data);

      newton_solver.verbose(m_verbose);

      // Precompute gradient if ok
      if(m_rp.m_precomputation){
         timecounter t1;
         t1.tic();
         this->pre_computation();
         t1.toc();
         if(m_verbose){
            std::cout << "-Precomputation: " << t1.to_double() << " sec" << std::endl;
         }
      }

      // List of time step
      std::list<time_step> list_step;

      scalar_type time1 = 0.0;
      for (size_t n = 0; n < m_rp.m_time_step.size(); n++)
      {
         const auto time_info = m_rp.m_time_step[n];
         const scalar_type time2 = time_info.first;
         const scalar_type delta_t = (time2 - time1)/time_info.second;
         for (size_t i = 0; i < time_info.second; i++) {
            time_step step;
            step.time = time1 + (i+1) * delta_t;
            step.level = 1;
            list_step.push_back(step);
         }
         time1 = time2;
      }

      size_t current_step = 0;
      size_t total_step = list_step.size();
      scalar_type old_time = 0.0;

      //time of saving
      bool time_saving(false);
      if(m_rp.m_n_time_save > 0.0){
         time_saving = true;
      }

      std::vector<std::array<scalar_type,2>> tab_traction;
      tab_traction.push_back({0,0});

      while(!list_step.empty())
      {
         current_step += 1;
         time_step step = list_step.front();
         const scalar_type current_time = step.time;
         if(m_verbose){
            std::cout << "--------------------------------------------------------------------------------------------------------------------------------" << std::endl;
            std::cout << "***************************************** Time : " << current_time << " sec (step: " << current_step
            << "/" << total_step << ", sublevel: " << step.level << " ) ****************************************|" << std::endl;
         }

         //Create a load function for each time step
         auto rlf = [&lf, &current_time](const Point& p) -> auto {
            return disk::mm_prod(current_time, lf(p));

         };

         auto rbcf = [&bcf, &current_time](const Point& p) -> auto {
            return disk::mm_prod(current_time, bcf(p));
         };


         auto rg = [&g, &current_time](const Point& p) -> auto {
            return disk::mm_prod(current_time, g(p));
         };

         // Compute the solution for the current step
         const NewtonSolverInfo newton_info = newton_solver.compute(rlf, rbcf, rg, m_gradient_precomputed);
         si.updateInfo(newton_info);

         if(m_verbose){
            std::cout << "** Time in this step " << newton_info.m_time_newton << " sec" << std::endl;
            std::cout << "**** Assembly time: " << newton_info.m_assembly_info.m_time_assembly << " sec" << std::endl;
            std::cout << "****** Gradient reconstruction: " << newton_info.m_assembly_info.m_time_gradrec << " sec" << std::endl;
            std::cout << "****** Stabilisation: " << newton_info.m_assembly_info.m_time_stab << " sec" << std::endl;
            std::cout << "****** Mechanical computation: " << newton_info.m_assembly_info.m_time_elem << " sec" << std::endl;
            std::cout << "       *** Behavior computation: " << newton_info.m_assembly_info.m_time_law << " sec" << std::endl;
            std::cout << "****** Static condensation: " << newton_info.m_assembly_info.m_time_statcond << " sec" << std::endl;
            std::cout << "****** Postprocess time: " << newton_info.m_assembly_info.m_time_postpro << " sec" << std::endl;
            std::cout << "**** Solver time: " << newton_info.m_solve_info.m_time_solve << " sec" << std::endl;
         }

         // Test Convergence
         m_convergence = newton_solver.test_convergence();

         // If no convergence
         if(!m_convergence){
            // The number of sublevel is reached
            if(step.level > m_rp.m_sublevel){
               std::cout << "***********************************************************" << std::endl;
               std::cout << "***** PROBLEM OF CONVERGENCE: We stop the calcul here *****" << std::endl;
               std::cout << "***********************************************************" << std::endl;
               break;
            }
            else
            {
               // The step is splitted in two
               if(m_verbose){
                  std::cout << "***********************************************************" << std::endl;
                  std::cout << "*****     NO CONVERGENCE: We split the time step     ******" << std::endl;
                  std::cout << "***********************************************************" << std::endl;
               }
               total_step += 1;
               current_step -= 1;
               time_step new_step;
               new_step.time = old_time + (current_time - old_time)/2.0;
               new_step.level = step.level + 1;
               list_step.push_front(new_step);
            }
         }
         else{
            // Delete the current step of the list
            old_time = current_time;
            list_step.pop_front();

            // Save the solution
            newton_solver.save_solutions(m_solution_cells, m_solution_faces, m_solution_lagr, m_solution_data);

            // Export the solution
            if(time_saving){
               if(m_rp.m_time_save.front()< old_time + 1E-5){
                  std::cout << "** Save results" << std::endl;
                  std::string name = "result_k" +  std::to_string(m_rp.m_cell_degree) + "_l" + std::to_string(m_rp.m_face_degree)
                  + "_g" + std::to_string(m_rp.m_grad_degree) + "_t" + std::to_string(old_time) + "_";
                  this->compute_discontinuous_displacement(name + "DEPL_disc.msh");
                  this->compute_continuous_displacement(name +"DEPL_cont.msh");
                  this->compute_deformed(name +"DEF.msh");
                  this->compute_J_GP(name +"J_GP.msh");
                  this->compute_continuous_J(name +"J_cont.msh");
                  this->compute_discontinuous_J(name +"J_disc.msh");
//                   try {
//                      this->compute_discontinuous_Prr(name +"Prr.msh", "Prr");
//                   }
//                   catch(const std::invalid_argument& ia){
//                      std::cerr << "Invalid argument: " << ia.what()  << " in Prr_disc" << std::endl;
//                   }
//                   try {
//                      this->compute_discontinuous_Prr(name +"Poo.msh", "Poo");
//                   }
//                   catch(const std::invalid_argument& ia){
//                      std::cerr << "Invalid argument: " << ia.what() << " in Prr_disc" << std::endl;
//                   }
                  try {
                     this->compute_discontinuous_VMIS(name +"VM_disc.msh");
                  }
                  catch(const std::invalid_argument& ia){
                     std::cerr << "Invalid argument: " << ia.what() << " in VMIS_disc" << std::endl;
                  }
                  try {
                     this->compute_continuous_VMIS(name +"VM_cont.msh");
                  }
                  catch(const std::invalid_argument& ia){
                     std::cerr << "Invalid argument: " << ia.what() << " in VMIS_cont" << std::endl;
                  }
                  try {
                     this->compute_VMIS_GP(name +"VM_GP.msh");
                  }
                  catch(const std::invalid_argument& ia){
                     std::cerr << "Invalid argument: " << ia.what() << " in VM_GP" << std::endl;
                  }

                  m_rp.m_time_save.pop_front();
                  if(m_rp.m_time_save.empty())  time_saving = false;
               }
            }
            //traction discrete
            //const auto traction = this->compute_traction_RT(1);
            tab_traction.push_back({current_time, this->compute_traction_Pk(3)});
         }
      }

      this->save_traction(tab_traction, "traction.dat");

      si.m_time_step = total_step;

      ttot.toc();
      si.m_time_solver = ttot.to_double();
      return si;
   }

   bool test_convergence() const {return m_convergence;}
   size_t getDofs() const {return total_dof_depl_static;}

   // PostProcessing

   // Compute different L2 error
   template<typename AnalyticalSolution>
   scalar_type
   compute_l2_error_displacement(const AnalyticalSolution& as) const
   {
      scalar_type err_dof = scalar_type{0.0};

      projector_type projk(m_bqd);

      size_t cell_i(0);

      for (auto& cl : m_msh)
      {
         const auto x = m_solution_cells.at(cell_i++);
         const dynamic_vector<scalar_type> true_dof = projk.compute_cell(m_msh, cl, as);
         const dynamic_vector<scalar_type> comp_dof = x.block(0,0,true_dof.size(), 1);
         const dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
         err_dof += diff_dof.dot(projk.cell_mm * diff_dof);
      }

      return sqrt(err_dof);
   }


   template<typename AnalyticalSolution>
   scalar_type
   compute_l2_error_gradient(const AnalyticalSolution& grad) const
   {
      scalar_type err_dof = scalar_type{0.0};

      projector_type projk(m_bqd);

      gradrec_type gradrec(m_bqd);

      size_t cell_i(0);

      for (auto& cl : m_msh)
      {
         const auto x = m_solution_data.at(cell_i);
         matrix_dynamic GT;
         if(m_rp.m_precomputation){
            GT = m_gradient_precomputed[cell_i];
         }
         else{
            gradrec.compute(m_msh, cl, false);
            GT = gradrec.oper;
         }
         const dynamic_vector<scalar_type> GTu = GT*x;
         const dynamic_vector<scalar_type> true_dof = projk.compute_cell_grad(m_msh, cl, grad);
         const dynamic_vector<scalar_type> comp_dof = GTu.block(0,0,true_dof.size(), 1);
         const dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
         err_dof += diff_dof.dot(projk.grad_mm * diff_dof);
         cell_i++;
      }

      return sqrt(err_dof);
   }


   template<typename AnalyticalSolution>
   scalar_type
   compute_l2_error_energy(const AnalyticalSolution& grad) const
   {

      gradrec_type gradrec(m_bqd);
      projector_type projk(m_bqd);

      const NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);

      size_t cell_i(0);
      scalar_type error_energy(0.0);

      for (auto& cl : m_msh)
      {
         const auto x = m_solution_data.at(cell_i);
         matrix_dynamic GT;
         if(m_rp.m_precomputation){
            GT = m_gradient_precomputed[cell_i];
         }
         else{
            gradrec.compute(m_msh, cl, false);
            GT = gradrec.oper;
         }
         const dynamic_vector<scalar_type> GTu = GT*x;
         const dynamic_vector<scalar_type> true_dof = projk.compute_cell_grad(m_msh, cl, grad);

         const auto grad_quadpoints = m_bqd.grad_quadrature.integrate(m_msh, cl);

         for (auto& qp : grad_quadpoints)
         {
            const auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, qp.point());
            const auto GT_iqn = disk::hho::eval_gradient(GTu, gphi);
            const auto FT_iqn = disk::mechanics::convertGtoF(GT_iqn);
            const scalar_type energy_comp = law.compute_energy(FT_iqn);

            const auto GT_true = disk::hho::eval_gradient(true_dof, gphi);
            const auto FT_true = disk::mechanics::convertGtoF(GT_true);
            const scalar_type energy_true = law.compute_energy(FT_true);

            error_energy += qp.weight() * std::pow(energy_true - energy_comp, 2.0);
         }
         cell_i++;
      }

      return sqrt(error_energy);
   }

   template<typename AnalyticalSolution>
   scalar_type
   compute_l2_error_PK1(const AnalyticalSolution& grad) const
   {
      gradrec_type gradrec(m_bqd);
      projector_type projk(m_bqd);

      const NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);

      size_t cell_i(0);
      scalar_type error_PK1(0.0);

      for (auto& cl : m_msh)
      {
         const auto x = m_solution_data.at(cell_i);
         matrix_dynamic GT;
         if(m_rp.m_precomputation){
            GT = m_gradient_precomputed[cell_i];
         }
         else{
            gradrec.compute(m_msh, cl, false);
            GT = gradrec.oper;
         }
         const dynamic_vector<scalar_type> GTu = GT*x;
         const dynamic_vector<scalar_type> true_dof = projk.compute_cell_grad(m_msh, cl, grad);

         const auto grad_quadpoints = m_bqd.grad_quadrature.integrate(m_msh, cl);

         for (auto& qp : grad_quadpoints)
         {
            const auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, qp.point());
            const auto GT_iqn = disk::hho::eval_gradient(GTu, gphi);
            const auto FT_iqn = disk::mechanics::convertGtoF(GT_iqn);
            const auto PK1_comp = law.compute_PK1(FT_iqn);

            const auto GT_true = disk::hho::eval_gradient(true_dof, gphi);
            const auto FT_true = disk::mechanics::convertGtoF(GT_true);
            const auto PK1_true = law.compute_PK1(FT_true);

            const auto PK1_diff = (PK1_true - PK1_comp).eval();

            error_PK1 += qp.weight() * disk::mm_prod(PK1_diff, PK1_diff);
         }
         cell_i++;
      }

      return sqrt(error_PK1);
   }

   template<typename AnalyticalSolution>
   void
   compute_l2_error_displacement_GP(const std::string& filename, const AnalyticalSolution& as) const
   {
      visu::Gmesh msh; //creta a mesh

      std::vector<visu::Data> data; //create data (not used)
      std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point
      size_t nb_node =  msh.getNumberofNodes();

      size_t cell_i(0);
      for (auto& cl : m_msh)
      {
         const vector_dynamic x = m_solution_cells.at(cell_i++);
         const auto qps = m_bqd.cell_quadrature.integrate(m_msh, cl);
         for (auto& qp : qps)
         {
            nb_node++;
            const auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, qp.point());

            static_vector<scalar_type, DIM> depl = static_vector<scalar_type, DIM>::Zero();

            // Compute displacement
            for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_bqd.m_cell_degree()).size(); i += DIM){
               for(size_t j = 0; j < DIM; j++){
                  depl(j) += phi.at(i+j)(j) * x(i+j);
               }
            }
            // True displacement
            const auto true_depl = as(qp.point()); // a voir et projeté

            // Create a node at gauss point
            const visu::Node snode = visu::convertPoint(qp.point(), nb_node);
            const std::vector<double> value(1, std::sqrt(disk::mm_prod(depl - true_depl, depl - true_depl))); // save the solution at gauss point

            // Save value
            const visu::SubData sdata(value, snode);
            subdata.push_back(sdata); // add subdata
         }
      }
      // Create and init a nodedata view
      visu::NodeData nodedata(1, 0.0, "error_depl_gp", data, subdata);
      // Save the view
      nodedata.saveNodeData(filename, msh);
   }


   void
   compute_discontinuous_displacement(const std::string& filename) const
   {
      const size_t cell_degree = m_bqd.cell_degree();
      const size_t num_cell_dofs = (m_bqd.cell_basis.range(0, cell_degree)).size();

      visu::Gmesh gmsh(DIM);
      auto storage = m_msh.backend_storage();

      std::vector<visu::Data> data; //create data (not used)
      const std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point

      size_t cell_i(0);
      size_t nb_nodes(0);
      for (auto& cl : m_msh)
      {
         const vector_dynamic x = m_solution_cells.at(cell_i++);
         const auto cell_nodes = visu::cell_nodes(m_msh, cl);
         std::vector<visu::Node> new_nodes;

         //loop on the nodes of the cell
         for (size_t i = 0; i < cell_nodes.size(); i++)
         {
            nb_nodes++;
            const auto point_ids = cell_nodes[i];
            const auto pt = storage->points[point_ids];

            const auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt);

            std::vector<scalar_type> depl(3, scalar_type{0});
            std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};

            // Add a node
            visu::init_coor(pt, coor);
            const visu::Node tmp_node(coor, nb_nodes, 0);
            new_nodes.push_back(tmp_node);
            gmsh.addNode(tmp_node);

            // Plot displacement at node
            for (size_t i = 0; i < num_cell_dofs; i += DIM){
               for(size_t j = 0; j < DIM; j++){
                  depl[j] += phi[i+j](j) * x(i+j); // a voir
               }
            }

            const visu::Data datatmp(nb_nodes, depl);
            data.push_back(datatmp);
         }
         // Add new element
         visu::add_element(gmsh, new_nodes);
      }

      // Create and init a nodedata view
      visu::NodeData nodedata(3, 0.0, "depl_node", data, subdata);

      // Save the view
      nodedata.saveNodeData(filename, gmsh);
   }


   void
   compute_continuous_displacement(const std::string& filename) const
   {
      const size_t cell_degree = m_bqd.cell_degree();
      const size_t num_cell_dofs = (m_bqd.cell_basis.range(0, cell_degree)).size();

      visu::Gmesh gmsh = visu::convertMesh(m_msh);
      auto storage = m_msh.backend_storage();

      const std::vector<double> vzero(3,0.0);
      const size_t nb_nodes(gmsh.getNumberofNodes());

      //first(number of data at this node), second(cumulated value)
      std::vector<std::pair<size_t, std::vector<scalar_type> > > value(nb_nodes, std::make_pair(0, vzero));

      size_t cell_i(0);
      for (auto& cl : m_msh)
      {
         vector_dynamic x = m_solution_cells.at(cell_i++);
         auto cell_nodes = visu::cell_nodes(m_msh, cl);

         // Loop on the nodes of the cell
         for (size_t i = 0; i < cell_nodes.size(); i++)
         {
            const auto point_ids = cell_nodes[i];
            const auto pt = storage->points[point_ids];

            const auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt);

            std::vector<scalar_type> depl(3, scalar_type{0});

            // Compute displacement at node
            for (size_t i = 0; i < num_cell_dofs; i += DIM){
               for(size_t j = 0; j < DIM; j++){
                  depl[j] += phi.at(i+j)(j) * x(i+j); // a voir
               }
            }

            // Add displacement at node
            value[point_ids].first++;
            for(size_t j=0; j < DIM; j++){
               value[point_ids].second[j] += depl[j];
            }
         }
      }

      std::vector<visu::Data> data; //create data
      std::vector<visu::SubData> subdata; //create subdata
      data.reserve(nb_nodes); // data has a size of nb_node

      // Compute the average value and save it
      for(size_t  i_node = 0; i_node < value.size(); i_node++){
         std::vector<double> tmp_value(3,0.0);
         for(size_t j=0; j < DIM; j++)
            tmp_value[j] = value[i_node].second[j]/ double(value[i_node].first);

         const visu::Data tmp_data(i_node + 1, tmp_value);
         data.push_back(tmp_data);
      }

      // Create and init a nodedata view
      visu::NodeData nodedata(3, 0.0, "depl_node", data, subdata);
      // Save the view
      nodedata.saveNodeData(filename, gmsh);
   }

   void
   compute_deformed_CONT(const std::string& filename) const
   {
      const size_t cell_degree = m_bqd.cell_degree();
      const size_t num_cell_dofs = (m_bqd.cell_basis.range(0, cell_degree)).size();

      visu::Gmesh gmsh(DIM);
      auto storage = m_msh.backend_storage();

      size_t cell_i(0);
      size_t nb_nodes(0);
      for (auto& cl : m_msh)
      {
         const vector_dynamic x = m_solution_cells.at(cell_i++);
         const auto cell_nodes = visu::cell_nodes(m_msh, cl);
         std::vector<visu::Node> new_nodes;

         // Loop on nodes of the cell
         for (size_t i = 0; i < cell_nodes.size(); i++)
         {
            nb_nodes++;
            const auto point_ids = cell_nodes[i];
            const auto pt = storage->points[point_ids];

            const auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt);

            std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};
            std::array<double, 3> depl = {double{0.0}, double{0.0}, double{0.0}};

            // Compute displacement
            visu::init_coor(pt, coor);
            for (size_t i = 0; i < num_cell_dofs; i += DIM){
               for(size_t j=0; j < DIM; j++){
                  depl[j] += phi.at(i+j)(j) * x(i+j);
               }
            }

            // Compute new coordinates
            for(size_t j=0; j < DIM; j++)
               coor[j] += depl[j];

            // Save node
            const visu::Node tmp_node(coor, nb_nodes, 0);
            new_nodes.push_back(tmp_node);
            gmsh.addNode(tmp_node);

         }
         // Add new element
         visu::add_element(gmsh, new_nodes);
      }
      // Save mesh
      gmsh.writeGmesh(filename, 2);
   }

   void
   compute_deformed(const std::string& filename) const
   {
      const size_t cell_degree = m_bqd.cell_degree();
      const size_t num_cell_dofs = (m_bqd.cell_basis.range(0, cell_degree)).size();

      visu::Gmesh gmsh(DIM);
      auto storage = m_msh.backend_storage();

      size_t cell_i(0);
      size_t nb_nodes(0);
      for (auto& cl : m_msh)
      {
         const vector_dynamic x = m_solution_cells.at(cell_i++);
         const auto cell_nodes = visu::cell_nodes(m_msh, cl);
         std::vector<visu::Node> new_nodes;

         // Loop on nodes of the cell
         for (size_t i = 0; i < cell_nodes.size(); i++)
         {
            nb_nodes++;
            const auto point_ids = cell_nodes[i];
            const auto pt = storage->points[point_ids];

            const auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt);

            std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};
            std::array<double, 3> depl = {double{0.0}, double{0.0}, double{0.0}};

            // Compute displacement
            visu::init_coor(pt, coor);
            for (size_t i = 0; i < num_cell_dofs; i += DIM){
               for(size_t j=0; j < DIM; j++){
                  depl[j] += phi.at(i+j)(j) * x(i+j);
               }
            }

            // Compute new coordinates
            for(size_t j=0; j < DIM; j++)
               coor[j] += depl[j];

            // Save node
            const visu::Node tmp_node(coor, nb_nodes, 0);
            new_nodes.push_back(tmp_node);
            gmsh.addNode(tmp_node);

         }
         // Add new element
         visu::add_element(gmsh, new_nodes);
      }
      // Save mesh
      gmsh.writeGmesh(filename, 2);
   }

   void
   compute_discontinuous_PK1(const std::string& filename) const
   {
      visu::Gmesh gmsh(DIM);
      auto storage = m_msh.backend_storage();

      gradrec_type gradrec(m_bqd);
      const NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);
      //const CavitationLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);

      std::vector<visu::Data> data; //create data (not used)
      std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point

      size_t cell_i(0);
      size_t nb_nodes(0);
      for (auto& cl : m_msh)
      {
         matrix_dynamic GT;
         if(m_rp.m_precomputation){
            GT = m_gradient_precomputed[cell_i];
         }
         else{
            gradrec.compute(m_msh, cl, false);
            GT = gradrec.oper;
         }
         const vector_dynamic GT_uTF = GT * m_solution_data.at(cell_i);

         const auto cell_nodes = visu::cell_nodes(m_msh, cl);
         std::vector<visu::Node> new_nodes;

         // Loop on nodes
         for (size_t i = 0; i < cell_nodes.size(); i++)
         {
            nb_nodes++;
            const auto point_ids = cell_nodes[i];
            const auto pt = storage->points[point_ids];

            const auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, pt);

            const auto GT_iqn = disk::hho::eval_gradient(GT_uTF, gphi);
            const auto FT_iqn = disk::mechanics::convertGtoF(GT_iqn);

            const auto PK1 = law.compute_PK1(FT_iqn);

            std::vector<scalar_type> PK1_pt(1, PK1.norm());
            std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};

            visu::init_coor(pt, coor);
            const visu::Node tmp_node(coor, nb_nodes, 0);
            new_nodes.push_back(tmp_node);
            gmsh.addNode(tmp_node);

            // Save data
            const visu::Data datatmp(nb_nodes, PK1_pt);
            data.push_back(datatmp);
         }
         // add new element
         visu::add_element(gmsh, new_nodes);
         cell_i++;
      }

      // Save
      visu::NodeData nodedata(1, 0.0, "PK1", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }

   // Von Mises stress postprocessing
   void
   compute_discontinuous_VMIS(const std::string& filename) const
   {
      visu::Gmesh gmsh(DIM);
      auto storage = m_msh.backend_storage();

      gradrec_type gradrec(m_bqd);
      //const NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);
      const CavitationLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);

      std::vector<visu::Data> data; //create data (not used)
      std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point

      size_t cell_i(0);
      size_t nb_nodes(0);
      for (auto& cl : m_msh)
      {
         matrix_dynamic GT;
         if(m_rp.m_precomputation){
            GT = m_gradient_precomputed[cell_i];
         }
         else{
            gradrec.compute(m_msh, cl, false);
            GT = gradrec.oper;
         }
         const vector_dynamic GT_uTF = GT * m_solution_data.at(cell_i);

         const auto cell_nodes = visu::cell_nodes(m_msh, cl);
         std::vector<visu::Node> new_nodes;

         // Loop on nodes
         for (size_t i = 0; i < cell_nodes.size(); i++)
         {
            nb_nodes++;
            const auto point_ids = cell_nodes[i];
            const auto pt = storage->points[point_ids];

            const auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, pt);

            const auto GT_iqn = disk::hho::eval_gradient(GT_uTF, gphi);
            const auto FT_iqn = disk::mechanics::convertGtoF(GT_iqn);

            const auto PK1 = law.compute_PK1(FT_iqn);

            const auto sigma = disk::mechanics::convertPK1toCauchy(PK1, FT_iqn);

            scalar_type vm(0.0);

            vm = sigma(0,0)*sigma(0,0) + sigma(1,1)*sigma(1,1) - sigma(0,0)*sigma(1,1) + 3*sigma(0,1)*sigma(0,1);

            if(DIM==3){
               vm += sigma(2,2)*sigma(2,2) - sigma(0,0)*sigma(2,2) - sigma(1,1)*sigma(2,2);
               vm += 3*(sigma(1,2)*sigma(1,2) + sigma(0,2)*sigma(0,2));
            }

            vm = sqrt(vm);

            std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};

            visu::init_coor(pt, coor);
            const visu::Node tmp_node(coor, nb_nodes, 0);
            new_nodes.push_back(tmp_node);
            gmsh.addNode(tmp_node);

            const std::vector<scalar_type> value(1,vm);
            const visu::Data datatmp(nb_nodes, value);
            data.push_back(datatmp);
         }
         // add new element
         visu::add_element(gmsh, new_nodes);
         cell_i++;
      }

      visu::NodeData nodedata(1, 0.0, "VM", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }

   void
   compute_continuous_VMIS(const std::string& filename) const
   {
      visu::Gmesh gmsh = visu::convertMesh(m_msh);
      auto storage = m_msh.backend_storage();

      const size_t nb_nodes(gmsh.getNumberofNodes());

      gradrec_type gradrec(m_bqd);

      //first(number of data at this node), second(cumulated value)
      std::vector<std::pair<size_t, scalar_type> > value(nb_nodes, std::make_pair(0, 0.0));

      //const NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);
      const CavitationLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);

      size_t cell_i(0);
      for (auto& cl : m_msh)
      {
         matrix_dynamic GT;
         if(m_rp.m_precomputation){
            GT = m_gradient_precomputed[cell_i];
         }
         else{
            gradrec.compute(m_msh, cl, false);
            GT = gradrec.oper;
         }
         const vector_dynamic GT_uTF = GT * m_solution_data.at(cell_i);
         const auto cell_nodes = visu::cell_nodes(m_msh, cl);

         // Loop on nodes
         for (size_t i = 0; i < cell_nodes.size(); i++)
         {
            const auto point_ids = cell_nodes[i];
            const auto pt = storage->points[point_ids];

            const auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, pt);

            const auto GT_iqn = disk::hho::eval_gradient(GT_uTF, gphi);
            const auto FT_iqn = disk::mechanics::convertGtoF(GT_iqn);

            const auto PK1 = law.compute_PK1(FT_iqn);

            const auto sigma = disk::mechanics::convertPK1toCauchy(PK1, FT_iqn);

            scalar_type vm(0.0);

            vm = sigma(0,0)*sigma(0,0) + sigma(1,1)*sigma(1,1) - sigma(0,0)*sigma(1,1) + 3*sigma(0,1)*sigma(0,1);

            if(DIM==3){
               vm += sigma(2,2)*sigma(2,2) - sigma(0,0)*sigma(2,2) - sigma(1,1)*sigma(2,2);
               vm += 3*(sigma(1,2)*sigma(1,2) + sigma(0,2)*sigma(0,2));
            }

            vm = sqrt(vm);

            // Add VM at node
            value[point_ids].first++;
            value[point_ids].second += vm;
         }

         cell_i++;
      }

      std::vector<visu::Data> data; //create data
      std::vector<visu::SubData> subdata; //create subdata
      data.reserve(nb_nodes); // data has a size of nb_node

      // Compute the average value and save it
      for(size_t  i_node = 0; i_node < value.size(); i_node++){
         const std::vector<double> tmp_value(1, value[i_node].second/ double(value[i_node].first));

         const visu::Data tmp_data(i_node + 1, tmp_value);
         data.push_back(tmp_data);
      }

      visu::NodeData nodedata(1, 0.0, "VM", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }


   void
   compute_VMIS_GP(const std::string& filename) const
   {
      visu::Gmesh gmsh = visu::convertMesh(m_msh);
      auto storage = m_msh.backend_storage();

      std::vector<visu::Data> data; //create data (not used)
      std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point
      size_t nb_nodes(gmsh.getNumberofNodes());

      gradrec_type gradrec(m_bqd);
      //const NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);
      const CavitationLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);

      size_t cell_i(0);
      for (auto& cl : m_msh)
      {
         matrix_dynamic GT;
         if(m_rp.m_precomputation){
            GT = m_gradient_precomputed[cell_i];
         }
         else{
            gradrec.compute(m_msh, cl, false);
            GT = gradrec.oper;
         }
         const vector_dynamic GT_uTF = GT * m_solution_data.at(cell_i);

         const auto qps = m_bqd.grad_quadrature.integrate(m_msh, cl);

         // Loop on nodes
         for (auto& qp : qps)
         {
            const auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, qp.point());
            const auto GT_iqn = disk::hho::eval_gradient(GT_uTF, gphi);
            const auto FT_iqn = disk::mechanics::convertGtoF(GT_iqn);

            const auto PK1 = law.compute_PK1(FT_iqn);

            const auto sigma = disk::mechanics::convertPK1toCauchy(PK1, FT_iqn);

            scalar_type vm(0.0);

            vm = sigma(0,0)*sigma(0,0) + sigma(1,1)*sigma(1,1) - sigma(0,0)*sigma(1,1) + 3*sigma(0,1)*sigma(0,1);

            if(DIM==3){
               vm += sigma(2,2)*sigma(2,2) - sigma(0,0)*sigma(2,2) - sigma(1,1)*sigma(2,2);
               vm += 3*(sigma(1,2)*sigma(1,2) + sigma(0,2)*sigma(0,2));
            }

            vm = sqrt(vm);

            // Add GP
            // Create a node at gauss point
            nb_nodes++;
            const visu::Node new_node = visu::convertPoint(qp.point(), nb_nodes);
            const std::vector<double> value(1, vm);
            const visu::SubData sdata(value, new_node);
            subdata.push_back(sdata); // add subdata
         }
         cell_i++;
      }

      // Save
      visu::NodeData nodedata(1, 0.0, "VM", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }


   template<typename AnalyticalSolution>
   void
   compute_l2error_VMIS_GP(const std::string& filename, const AnalyticalSolution& grad) const
   {
      visu::Gmesh gmsh = visu::convertMesh(m_msh);
      auto storage = m_msh.backend_storage();

      std::vector<visu::Data> data; //create data (not used)
      std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point
      size_t nb_nodes(gmsh.getNumberofNodes());

      gradrec_type gradrec(m_bqd);
      const NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);
      //const CavitationLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);

      size_t cell_i(0);
      for (auto& cl : m_msh)
      {
         matrix_dynamic GT;
         if(m_rp.m_precomputation){
            GT = m_gradient_precomputed[cell_i];
         }
         else{
            gradrec.compute(m_msh, cl, false);
            GT = gradrec.oper;
         }
         const vector_dynamic GT_uTF = GT * m_solution_data.at(cell_i);

         const auto qps = m_bqd.grad_quadrature.integrate(m_msh, cl);

         // Loop on nodes
         for (auto& qp : qps)
         {
            const auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, qp.point());
            const auto GT_iqn = disk::hho::eval_gradient(GT_uTF, gphi);
            const auto GT_true = grad(qp.point());
            const auto FT_iqn = disk::mechanics::convertGtoF(GT_iqn);
            const auto FT_true = disk::mechanics::convertGtoF(GT_true);

            const auto PK1 = law.compute_PK1(FT_iqn);
            const auto PK1_true = law.compute_PK1(FT_true);

            const auto sigma = disk::mechanics::convertPK1toCauchy(PK1, FT_iqn);
            const auto sigma_true = disk::mechanics::convertPK1toCauchy(PK1_true, FT_true);

            scalar_type vm(0.0);

            vm = sigma(0,0)*sigma(0,0) + sigma(1,1)*sigma(1,1) - sigma(0,0)*sigma(1,1) + 3*sigma(0,1)*sigma(0,1);

            if(DIM==3){
               vm += sigma(2,2)*sigma(2,2) - sigma(0,0)*sigma(2,2) - sigma(1,1)*sigma(2,2);
               vm += 3*(sigma(1,2)*sigma(1,2) + sigma(0,2)*sigma(0,2));
            }

            vm = sqrt(vm);

            scalar_type vm_true(0.0);

            vm_true = sigma_true(0,0)*sigma_true(0,0) + sigma_true(1,1)*sigma_true(1,1)
            - sigma_true(0,0)*sigma_true(1,1) + 3*sigma_true(0,1)*sigma_true(0,1);

            if(DIM==3){
               vm_true += sigma_true(2,2)*sigma_true(2,2) - sigma_true(0,0)*sigma_true(2,2) - sigma_true(1,1)*sigma_true(2,2);
               vm_true += 3*(sigma_true(1,2)*sigma_true(1,2) + sigma_true(0,2)*sigma_true(0,2));
            }

            vm_true = sqrt(vm_true);

            // Add GP
            // Create a node at gauss point
            nb_nodes++;
            const visu::Node new_node = visu::convertPoint(qp.point(), nb_nodes);
            const std::vector<double> value(1, std::abs(vm - vm_true));
            const visu::SubData sdata(value, new_node);
            subdata.push_back(sdata); // add subdata
         }
         cell_i++;
      }

      // Save
      visu::NodeData nodedata(1, 0.0, "VM", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }

   template<typename AnalyticalSolution>
   void
   plot_analytical_VMIS(const std::string& filename, const AnalyticalSolution& grad) const
   {
      visu::Gmesh gmsh = visu::convertMesh(m_msh);
      auto storage = m_msh.backend_storage();

      const size_t nb_nodes(gmsh.getNumberofNodes());
      std::vector<visu::Data> data; //create data
      std::vector<visu::SubData> subdata; //create subdata
      data.reserve(nb_nodes); // data has a size of nb_node

      //first(number of data at this node), second(cumulated value)
      std::vector<std::pair<size_t, scalar_type> > value(nb_nodes, std::make_pair(0, 0.0));

      const NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);
      //const CavitationLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);

      size_t i_node(0);
      // Loop on nodes
      for (auto& pt : storage->points)
      {
         const auto GT = grad(pt);
         const auto FT = disk::mechanics::convertGtoF(GT);

         const auto PK1 = law.compute_PK1(FT);

         const auto sigma = disk::mechanics::convertPK1toCauchy(PK1, FT);

         scalar_type vm(0.0);

         vm = sigma(0,0)*sigma(0,0) + sigma(1,1)*sigma(1,1) - sigma(0,0)*sigma(1,1) + 3*sigma(0,1)*sigma(0,1);

         if(DIM==3){
            vm += sigma(2,2)*sigma(2,2) - sigma(0,0)*sigma(2,2) - sigma(1,1)*sigma(2,2);
            vm += 3*(sigma(1,2)*sigma(1,2) + sigma(0,2)*sigma(0,2));
         }

         vm = sqrt(vm);

         // Add VM at node
         const std::vector<double> tmp_value(1, double(vm));
         const visu::Data tmp_data(i_node + 1, tmp_value);
         data.push_back(tmp_data);
         i_node++;
      }

      visu::NodeData nodedata(1, 0.0, "VM", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }

   // Jacobian J PostProcessing
   void
   compute_discontinuous_J(const std::string& filename) const
   {
      visu::Gmesh gmsh(DIM);
      auto storage = m_msh.backend_storage();

      gradrec_type gradrec(m_bqd);
      //const NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);
      const CavitationLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);

      std::vector<visu::Data> data; //create data (not used)
      std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point

      size_t cell_i(0);
      size_t nb_nodes(0);
      for (auto& cl : m_msh)
      {
         matrix_dynamic GT;
         if(m_rp.m_precomputation){
            GT = m_gradient_precomputed[cell_i];
         }
         else{
            gradrec.compute(m_msh, cl, false);
            GT = gradrec.oper;
         }
         const vector_dynamic GT_uTF = GT * m_solution_data.at(cell_i);

         const auto cell_nodes = visu::cell_nodes(m_msh, cl);
         std::vector<visu::Node> new_nodes;

         // Loop on nodes
         for (size_t i = 0; i < cell_nodes.size(); i++)
         {
            nb_nodes++;
            const auto point_ids = cell_nodes[i];
            const auto pt = storage->points[point_ids];

            const auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, pt);

            const auto GT_iqn = disk::hho::eval_gradient(GT_uTF, gphi);
            const auto FT_iqn = disk::mechanics::convertGtoF(GT_iqn);
            const auto J_iqn = FT_iqn.determinant();

            std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};

            visu::init_coor(pt, coor);
            const visu::Node tmp_node(coor, nb_nodes, 0);
            new_nodes.push_back(tmp_node);
            gmsh.addNode(tmp_node);

            const std::vector<scalar_type> value(1, J_iqn);
            const visu::Data datatmp(nb_nodes, value);
            data.push_back(datatmp);
         }
         // add new element
         visu::add_element(gmsh, new_nodes);
         cell_i++;
      }

      visu::NodeData nodedata(1, 0.0, "J", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }

   void
   compute_continuous_J(const std::string& filename) const
   {
      visu::Gmesh gmsh = visu::convertMesh(m_msh);
      auto storage = m_msh.backend_storage();

      const size_t nb_nodes(gmsh.getNumberofNodes());

      gradrec_type gradrec(m_bqd);

      //first(number of data at this node), second(cumulated value)
      std::vector<std::pair<size_t, scalar_type> > value(nb_nodes, std::make_pair(0, 0.0));

      //const NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);
      const CavitationLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);

      size_t cell_i(0);
      for (auto& cl : m_msh)
      {
         matrix_dynamic GT;
         if(m_rp.m_precomputation){
            GT = m_gradient_precomputed[cell_i];
         }
         else{
            gradrec.compute(m_msh, cl, false);
            GT = gradrec.oper;
         }
         const vector_dynamic GT_uTF = GT * m_solution_data.at(cell_i);
         const auto cell_nodes = visu::cell_nodes(m_msh, cl);

         // Loop on nodes
         for (size_t i = 0; i < cell_nodes.size(); i++)
         {
            const auto point_ids = cell_nodes[i];
            const auto pt = storage->points[point_ids];

            const auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, pt);

            const auto GT_iqn = disk::hho::eval_gradient(GT_uTF, gphi);
            const auto FT_iqn = disk::mechanics::convertGtoF(GT_iqn);
            const auto J_iqn = FT_iqn.determinant();

            // Add VM at node
            value[point_ids].first++;
            value[point_ids].second += J_iqn;
         }

         cell_i++;
      }

      std::vector<visu::Data> data; //create data
      std::vector<visu::SubData> subdata; //create subdata
      data.reserve(nb_nodes); // data has a size of nb_node

      // Compute the average value and save it
      for(size_t  i_node = 0; i_node < value.size(); i_node++){
         const std::vector<double> tmp_value(1, value[i_node].second/ double(value[i_node].first));

         const visu::Data tmp_data(i_node + 1, tmp_value);
         data.push_back(tmp_data);
      }

      visu::NodeData nodedata(1, 0.0, "J", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }


   void
   compute_J_GP(const std::string& filename) const
   {
      visu::Gmesh gmsh = visu::convertMesh(m_msh);
      auto storage = m_msh.backend_storage();

      std::vector<visu::Data> data; //create data (not used)
      std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point
      size_t nb_nodes(gmsh.getNumberofNodes());

      gradrec_type gradrec(m_bqd);
      //const NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);
      const CavitationLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);

      size_t cell_i(0);
      for (auto& cl : m_msh)
      {
         matrix_dynamic GT;
         if(m_rp.m_precomputation){
            GT = m_gradient_precomputed[cell_i];
         }
         else{
            gradrec.compute(m_msh, cl, false);
            GT = gradrec.oper;
         }
         const vector_dynamic GT_uTF = GT * m_solution_data.at(cell_i);

         const auto qps = m_bqd.grad_quadrature.integrate(m_msh, cl);

         // Loop on nodes
         for (auto& qp : qps)
         {
            const auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, qp.point());
            const auto GT_iqn = disk::hho::eval_gradient(GT_uTF, gphi);
            const auto FT_iqn = disk::mechanics::convertGtoF(GT_iqn);
            const auto J_iqn = FT_iqn.determinant();

            // Add GP
            // Create a node at gauss point
            nb_nodes++;
            const visu::Node new_node = visu::convertPoint(qp.point(), nb_nodes);
            const std::vector<double> value(1, J_iqn);
            const visu::SubData sdata(value, new_node);
            subdata.push_back(sdata); // add subdata
         }
         cell_i++;
      }

      // Save
      visu::NodeData nodedata(1, 0.0, "J_GP", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }


   template<typename AnalyticalSolution>
   void
   compute_l2error_J_GP(const std::string& filename, const AnalyticalSolution& grad) const
   {
      visu::Gmesh gmsh = visu::convertMesh(m_msh);
      auto storage = m_msh.backend_storage();

      std::vector<visu::Data> data; //create data (not used)
      std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point
      size_t nb_nodes(gmsh.getNumberofNodes());

      gradrec_type gradrec(m_bqd);
      const NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);
      //const CavitationLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);

      size_t cell_i(0);
      for (auto& cl : m_msh)
      {
         matrix_dynamic GT;
         if(m_rp.m_precomputation){
            GT = m_gradient_precomputed[cell_i];
         }
         else{
            gradrec.compute(m_msh, cl, false);
            GT = gradrec.oper;
         }
         const vector_dynamic GT_uTF = GT * m_solution_data.at(cell_i);

         const auto qps = m_bqd.grad_quadrature.integrate(m_msh, cl);

         // Loop on nodes
         for (auto& qp : qps)
         {
            const auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, qp.point());
            const auto GT_iqn = disk::hho::eval_gradient(GT_uTF, gphi);
            const auto GT_true = grad(qp.point());
            const auto FT_iqn = disk::mechanics::convertGtoF(GT_iqn);
            const auto FT_true = disk::mechanics::convertGtoF(GT_true);
            const auto J_iqn = FT_iqn.determinant();
            const auto J_true = FT_true.determinant();

            // Add GP
            // Create a node at gauss point
            nb_nodes++;
            const visu::Node new_node = visu::convertPoint(qp.point(), nb_nodes);
            const std::vector<double> value(1, std::abs(J_iqn - J_true));
            const visu::SubData sdata(value, new_node);
            subdata.push_back(sdata); // add subdata
         }
         cell_i++;
      }

      // Save
      visu::NodeData nodedata(1, 0.0, "J_GP", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }

   template<typename AnalyticalSolution>
   void
   plot_analytical_J(const std::string& filename, const AnalyticalSolution& grad) const
   {
      visu::Gmesh gmsh = visu::convertMesh(m_msh);
      auto storage = m_msh.backend_storage();

      const size_t nb_nodes(gmsh.getNumberofNodes());
      std::vector<visu::Data> data; //create data
      std::vector<visu::SubData> subdata; //create subdata
      data.reserve(nb_nodes); // data has a size of nb_node

      //first(number of data at this node), second(cumulated value)
      std::vector<std::pair<size_t, scalar_type> > value(nb_nodes, std::make_pair(0, 0.0));

      const NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);
      //const CavitationLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);

      size_t i_node(0);
      // Loop on nodes
      for (auto& pt : storage->points)
      {
         const auto GT = grad(pt);
         const auto FT = convertGtoF(GT);
         const auto J = FT.determinant();

         // Add VM at node
         const std::vector<double> tmp_value(1, double(J));
         const visu::Data tmp_data(i_node + 1, tmp_value);
         data.push_back(tmp_data);
         i_node++;
      }

      visu::NodeData nodedata(1, 0.0, "J", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }


   /// PROVISOIRE A SUPPRIMER
   // A supprimer
   //compute PK in cylindrical base
   void
   compute_discontinuous_Prr(const std::string& filename, const std::string compo = "Prr") const
   {
      visu::Gmesh gmsh(DIM);
      auto storage = m_msh.backend_storage();

      gradrec_type gradrec(m_bqd);
      const NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);
      //const CavitationLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);

      std::vector<visu::Data> data; //create data (not used)
      std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point

      size_t cell_i(0);
      size_t nb_nodes(0);
      for (auto& cl : m_msh)
      {
         gradrec.compute(m_msh, cl, false);
         const vector_dynamic GT_uTF = gradrec.oper * m_solution_data.at(cell_i);

         auto cell_nodes = visu::cell_nodes(m_msh, cl);
         std::vector<visu::Node> new_nodes;
         for (size_t i = 0; i < cell_nodes.size(); i++)
         {
            nb_nodes++;
            auto point_ids = cell_nodes[i];
            auto pt = storage->points[point_ids];

            auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, pt);

            const auto GT_iqn = disk::hho::eval_gradient(GT_uTF, gphi);
            const auto FT_iqn = disk::mechanics::convertGtoF(GT_iqn);

            auto PK1= law.compute_PK1(FT_iqn);

            vector_dynamic er;
            er.resize(DIM);
            er(0) = pt.x(); er(1) = pt.y();
            if(DIM==3) er(2) = 0.0;
            er /= er.norm();
            vector_dynamic eO;
            eO.resize(DIM);
            eO(0) = -er(1); eO(1) = er(0);
            if(DIM==3) eO(2) = 0.0;

            std::vector<scalar_type> PK1rr(1,0.0) ;
            if(compo == "Prr")
               PK1rr[0] = er.dot(PK1*er);
            else if(compo == "Poo")
               PK1rr[0] = eO.dot(PK1*eO);
            else if(compo == "Pro")
               PK1rr[0] = er.dot(PK1*eO);
            else if(compo == "Por")
               PK1rr[0] = eO.dot(PK1*er);

            std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};

            visu::init_coor(pt, coor);
            visu::Node tmp_node(coor, nb_nodes, 0);
            new_nodes.push_back(tmp_node);
            gmsh.addNode(tmp_node);


            visu::Data datatmp(nb_nodes, PK1rr);
            data.push_back(datatmp);
         }
         // add new element
         visu::add_element(gmsh, new_nodes);
         cell_i++;
      }

      visu::NodeData nodedata(1, 0.0, "PK1", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }

   scalar_type
   compute_traction_RT(const size_t& id) const
   {
      Hyperelasticity::discrete_traction_bq<bqdata_type> traction(m_bqd);
      size_t cell_i = 0;

      scalar_type ret(0);
      scalar_type mesure(0);

      for (auto& cl : m_msh)
      {
         const auto fcs = faces(m_msh, cl);

         for(auto fc : fcs)
         {
            if(m_msh.is_boundary(fc)){
               const auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), fc);
               if (!eid.first)
                  throw std::invalid_argument("This is a bug: face not found");

               const auto face_id = eid.second;
               const size_t b_id = m_msh.boundary_id(face_id);

               if(b_id == id){
                  traction.proj_grad_space(m_msh, cl, m_gradient_precomputed[cell_i],
                                 m_solution_data[cell_i], m_elas_param);

                  //GP loop
                  const auto face_quadpoints = m_bqd.face_quadrature.integrate(m_msh, fc);

                  for (auto& qp : face_quadpoints)
                  {
                     const auto T_iqn = traction.eval_RT(m_msh, cl, fc, qp.point());
                     const auto n = disk::normal(m_msh, cl, fc);
                     ret += qp.weight() * std::pow(disk::mm_prod(T_iqn, n),2.0);
                  }

                  mesure += disk::measure(m_msh, fc);
               }
            }
         }

         cell_i++;
      }
      return std::sqrt(ret);
   }

   scalar_type
   compute_traction_Pk(const size_t& id) const
   {
      Hyperelasticity::discrete_traction_bq<bqdata_type> traction(m_bqd);
      size_t cell_i = 0;

      scalar_type ret(0);
      scalar_type mesure(0);

      for (auto& cl : m_msh)
      {
         const auto fcs = faces(m_msh, cl);

         for(size_t face_i = 0; face_i < fcs.size(); face_i++)
         {
            const auto fc = fcs[face_i];
            if(m_msh.is_boundary(fc)){
               const auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), fc);
               if (!eid.first)
                  throw std::invalid_argument("This is a bug: face not found");

               const auto face_id = eid.second;
               const size_t b_id = m_msh.boundary_id(face_id);

               if(b_id == id){
                  traction.proj_grad_space(m_msh, cl, m_gradient_precomputed[cell_i],
                                   m_solution_data[cell_i], m_elas_param);

                  if(m_bqd.cell_degree() == m_bqd.grad_degree()){
                     traction.compute_Sadjoint(m_msh, cl, m_solution_data[cell_i]);
                  }
                  //GP loop
                  const auto face_quadpoints = m_bqd.face_quadrature.integrate(m_msh, fc);

                  for (auto& qp : face_quadpoints)
                  {
                     if(m_bqd.cell_degree() == m_bqd.grad_degree()){
                        const auto T_iqn = traction.eval_Pk_unstable(m_msh, cl, face_i, qp.point(), m_rp.m_beta);
                        const auto n = disk::normal(m_msh, cl, fc);
                        ret += qp.weight() * std::abs(disk::mm_prod(T_iqn, n));
                     }
                     else{
                        const auto T_iqn = traction.eval_Pk_stable(m_msh, cl, fc, qp.point());
                        const auto n = disk::normal(m_msh, cl, fc);
                        ret += qp.weight() * std::abs(disk::mm_prod(T_iqn, n));
                     }
                  }

                  mesure += disk::measure(m_msh, fc);
               }
            }
         }

         cell_i++;
      }
      return std::sqrt(ret);
   }

   void save_traction(const std::vector<std::array<scalar_type,2>>& resu, const std::string& filename)
   {
      std::ofstream output;
      output.open(filename, std::ofstream::out | std::ofstream::app);

      if (!output.is_open())
      {
         std::cerr << "Unable to open file " << filename << std::endl;
      }

      output << "time" << "\t" << "traction" << std::endl;

      for(size_t i = 0; i < resu.size(); i++)
      {
         output << "(" << (resu[i])[0] << " , " << (resu[i])[1] << ")";
      }


      output.close();
   }


   // A supprimer
   std::pair<scalar_type,scalar_type>
   compute_l2_error_annulus(const std::string& file_error) const
   {
      std::ifstream   ifs(file_error);
      std::string     keyword;

      if (!ifs.is_open())
      {
         std::cout << "Error opening " << file_error << std::endl;
      }

      //ne sert a rien
      ifs >> keyword >> keyword >> keyword >> keyword >> keyword;
      size_t num(0);
      ifs >> num;

      matrix_dynamic mat;
      mat.resize(5, num);

      //in the order R, phi, dphi, Prr, Poo
      for (size_t i = 0; i < num; i++) {
         ifs >> mat(0,i) >> mat(1,i) >> mat(2,i) >> mat(3,i) >> mat(4,i);
      }

      ifs.close();

      gradrec_type gradrec(m_bqd);

      size_t i = 0;
      scalar_type error_depl = 0.0;
      scalar_type error_grad = 0.0;

      for (auto& cl : m_msh)
      {
         auto x = m_solution_data.at(i++);
         gradrec.compute(m_msh, cl, false);
         dynamic_vector<scalar_type> GTu = gradrec.oper*x;

         auto grad_quadpoints = m_bqd.grad_quadrature.integrate(m_msh, cl);

         for (auto& qp : grad_quadpoints)
         {
            //compute depl
            vector_dynamic depl; depl.resize(2); depl.setConstant(0.0);
            auto c_phi = m_bqd.cell_basis.eval_functions(m_msh, cl, qp.point());
            for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_bqd.cell_degree()).size(); i += DIM)
               for(size_t j=0; j < DIM; j++)
                  depl[j] += c_phi.at(i+j)(j) * x(i+j); // a voir

                  // compute grad_quadpoints
            const auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, qp.point());
            const auto GT_iqn = disk::hho::eval_gradient(GTu, gphi);
            const auto FT_iqn = disk::mechanics::convertGtoF(GT_iqn);
            //compute er and eo
            vector_dynamic er; er.resize(2); er(0) = qp.point().x(); er(1) = qp.point().y();
            const scalar_type R = er.norm();
            er /= R;
            vector_dynamic eo; eo.resize(2); eo(0) = -er(1); eo(1) = er(0);

            const scalar_type dphi = er.dot(FT_iqn*er);
            const scalar_type phi = R*eo.dot(FT_iqn*eo);
            const scalar_type cphi = depl.dot(er);
            //
            size_t ind_R(0);
            for (size_t i = 0; i < num; i++) {
               if(R < mat(0,i)){
                  ind_R = i;
                  break;
               }
               if(i == num-1){
                  ind_R = i;
               }
            }

            // interpolation lineaire du deplacement
            const scalar_type a = (mat(1,ind_R) - mat(1,ind_R -1)) / (mat(0,ind_R) - mat(0,ind_R -1)) ;
            const scalar_type b = mat(1,ind_R -1) - a * mat(0,ind_R -1);
            const scalar_type phi_R = a*R + b;
            const scalar_type depl_R = phi_R -R;

            // interpolation lineaire du gradient

            const scalar_type c = (mat(2,ind_R) - mat(2,ind_R -1)) / (mat(0,ind_R) - mat(0,ind_R -1)) ;
            const scalar_type d = mat(2,ind_R -1) - c * mat(0,ind_R -1);
            const scalar_type dphi_R = c*R + d;
            matrix_dynamic grad_ref; grad_ref.resize(2,2); grad_ref.setConstant(0.0);

            //compute L2 error depl

            const scalar_type relative_displ = (depl_R- cphi);

            error_depl += qp.weight() * relative_displ * relative_displ;

            //compute l2 error gradient
            error_grad += qp.weight() * (std::pow(dphi_R - dphi, 2.0) + std::pow(phi_R - phi, 2.0));
         }
      }

      error_depl = sqrt(error_depl);
      error_grad = sqrt(error_grad);

      std::cout << "ERROR L2 ANNULUS" << '\n';
      std::cout << "L2 DEPL: " << error_depl << '\n';
      std::cout << "L2 GRAD: " << error_grad << '\n';

      return std::make_pair(error_depl,error_grad);
   }

   // A supprimer
   std::array<scalar_type, 3>
   displacement_node(const size_t num_node) const
   {
      auto storage = m_msh.backend_storage();
      std::array<scalar_type, 3> depl = {double{0.0}, double{0.0}, double{0.0}};
      //avooir cell cell
      size_t cell_i = 0;
      for (auto& cl : m_msh)
      {
         vector_dynamic x = m_solution_cells.at(cell_i++);
         auto cell_nodes = visu::cell_nodes(m_msh, cl);
         for (size_t i = 0; i < cell_nodes.size(); i++)
         {
            auto point_ids = cell_nodes[i];
            if(point_ids == num_node){
               auto pt = storage->points[point_ids];

               auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt);

               // plot magnitude at node
               for (size_t k = 0; k < m_bqd.cell_basis.range(0, m_bqd.cell_degree()).size(); k += DIM)
                  for(size_t j=0; j < DIM; j++)
                     depl[j] += phi.at(k+j)(j) * x(k+j); // a voir

                     return depl;
            }
         }
      }

      std::cout << "Invalid node number" << std::endl;
      return depl;
   }

};
