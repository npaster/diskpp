/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016, 2017 - matteo.cicuttin@enpc.fr
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

#include <iostream>
#include <algorithm>
#include <vector>

#include <sstream>

#include "../../config.h"

#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "hho/hho.hpp"
#include "hho/hho_bq.hpp"
#include "hho/hho_vector.hpp"
#include "bases/bases_utils.hpp"
#include "NewtonSolver/newton_solver.hpp"
#include "ElasticityParameters.hpp"
#include "BoundaryConditions.hpp"
#include "Informations.hpp"
#include "Parameters.hpp"

#include "../exemple_visualisation/visualisation/gmshDisk.hpp"
#include "../exemple_visualisation/visualisation/gmshConvertMesh.hpp"

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

   typedef disk::basis_quadrature_data_full<mesh_type, disk::scaled_monomial_vector_basis,
                                                      disk::scaled_monomial_matrix_basis,
                                                      disk::quadrature> bqdata_type;

   typedef dynamic_matrix<scalar_type>         matrix_dynamic;
   typedef dynamic_vector<scalar_type>         vector_dynamic;

   typedef disk::gradient_reconstruction_elas_full_bq<bqdata_type>     gradrec_type;
   typedef disk::projector_elas_bq<bqdata_type>                        projector_type;

   bqdata_type     m_bqd;

   const mesh_type& m_msh;

   std::vector<vector_dynamic>         m_solution_data;
   std::vector<vector_dynamic>         m_solution_cells, m_solution_faces, m_solution_lagr;

   bool m_verbose, m_convergence;

   ElasticityParameters m_elas_param;
   param_type m_rp;

   BoundaryConditions m_boundary_condition;

   size_t total_dof_depl_static;

public:
   hyperelasticity_solver(const mesh_type& msh, const param_type& rp, const ElasticityParameters elas_param,
                         const std::vector<BoundaryType>& boundary_neumann,
                         const std::vector<BoundaryType>& boundary_dirichlet)
   : m_msh(msh), m_verbose(rp.m_verbose), m_convergence(false), m_elas_param(elas_param), m_rp(rp)
   {
      size_t face_degree = rp.m_face_degree;
      if(face_degree <= 0)
      {
         std::cout << "'face_degree' should be > 0. Reverting to 1." << std::endl;
         face_degree = 1;
      }

      m_rp.m_face_degree = face_degree;

      size_t cell_degree = rp.m_face_degree + rp.m_l;
      if(face_degree-1  > cell_degree or cell_degree > face_degree +1 )
      {
         std::cout << "'cell_degree' should be 'face_degree + 1' => 'cell_degree' => 'face_degree -1'. Reverting to 'face_degree'." << std::endl;
         cell_degree = face_degree;
      }

      m_rp.m_cell_degree = cell_degree;

      size_t grad_degree = rp.m_grad_degree;
      if(grad_degree  < cell_degree)
      {
         std::cout << "'grad_degree' should be > 'cell_degree'. Reverting to 'cell_degree'." << std::endl;
         grad_degree = cell_degree;
      }

      m_rp.m_grad_degree = grad_degree;

      m_bqd = bqdata_type(face_degree, cell_degree, grad_degree);

      if(m_verbose){
         m_bqd.info_degree();
         m_rp.infos();
      }

      m_boundary_condition = BoundaryConditions(msh, boundary_neumann, boundary_dirichlet);

   }

   bool    verbose(void) const     { return m_verbose; }
   void    verbose(bool v)         { m_verbose = v; }



   //template<typename DeplFunction, typename StressFunction>
   void
   compute_initial_state()//const DeplFunction& df, const StressFunction& bcf)
   {
      m_solution_data.clear();
      m_solution_cells.clear();
      m_solution_faces.clear();
      m_solution_lagr.clear();

      const size_t nb_faces_dirichlet = m_boundary_condition.nb_faces_dirichlet();
      const size_t nb_lag_conditions = m_boundary_condition.nb_lag();

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
         auto fcs = faces(m_msh, cl);
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


      if(m_verbose){
         std::cout << "** Numbers of cells: " << m_msh.cells_size()  << std::endl;
         std::cout << "** Numbers of faces: " << m_msh.faces_size() << " ( boundary faces: "
         <<  m_msh.boundary_faces_size() << " )"  << std::endl;
         std::cout << "** Numbers of dofs: " << total_dof + total_lagr  << std::endl;
         std::cout << "** including " << total_lagr << " Lagrange multipliers"  << std::endl;
         std::cout << "** After static condensation: "  << std::endl;
         std::cout << "** Numbers of dofs: " << m_msh.faces_size() * num_face_dofs + total_lagr  << std::endl;
         std::cout << "** including " << total_lagr << " Lagrange multipliers"  << std::endl;
      }

      total_dof_depl_static = m_msh.faces_size() * num_face_dofs;
      //provisoire

//       for(size_t i = 0; i < m_msh.cells_size(); i++)
//       {
//          m_solution_data[i].setConstant(10.0);
//          m_solution_cells[i].setConstant(10.0);
//       }
//
//       for(size_t i = 0; i < m_msh.faces_size(); i++)
//       {
//          m_solution_faces[i].setConstant(10.0);
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

      NewtonRaphson_solver_hyperelasticity<bqdata_type>
         newton_solver(m_msh, m_bqd, m_rp, m_elas_param, m_boundary_condition);

      newton_solver.initialize(m_solution_cells, m_solution_faces,
                                 m_solution_lagr, m_solution_data);

      newton_solver.verbose(m_verbose);

      std::list<time_step> list_step;

      scalar_type time1 = 0.0;
      scalar_type time2 = 0.0;
      for (size_t n = 0; n < m_rp.m_time_step.size(); n++)
      {
         auto time_info = m_rp.m_time_step[n];
         time2 = time_info.first;
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

      while(!list_step.empty())
      {
         current_step += 1;
         time_step step = list_step.front();
         const scalar_type current_time = step.time;
         if(m_verbose){
            std::cout << "--------------------------------------------------------------" << std::endl;
            std::cout << "*************** Time : " << current_time << " sec (step: " << current_step
            << "/" << total_step << ") *****************|" << std::endl;
         }

         auto rlf = [&lf, &current_time](const Point& p) -> auto {
            return disk::mm_prod(current_time, lf(p));

         };

         auto rbcf = [&bcf, &current_time](const Point& p) -> auto {
            return disk::mm_prod(current_time, bcf(p));
         };


         auto rg = [&g, &current_time](const Point& p) -> auto {
            return disk::mm_prod(current_time, g(p));
         };

         NewtonSolverInfo newton_info = newton_solver.compute(rlf, rbcf, rg);
         si.updateInfo(newton_info);

         if(m_verbose){
            std::cout << "** Time in this step " << newton_info.m_time_newton << " sec" << std::endl;
            std::cout << "**** Assembly time: " << newton_info.m_assembly_info.m_time_assembly << " sec" << std::endl;
            std::cout << "****** Gradient reconstruction: " << newton_info.m_assembly_info.m_time_gradrec << " sec" << std::endl;
            std::cout << "****** Stabilisation: " << newton_info.m_assembly_info.m_time_stab << " sec" << std::endl;
            std::cout << "****** Elementary computation: " << newton_info.m_assembly_info.m_time_elem << " sec" << std::endl;
            std::cout << "       *** Behavior computation: " << newton_info.m_assembly_info.m_time_law << " sec" << std::endl;
            std::cout << "       *** Adaptative stabilization: " << newton_info.m_assembly_info.m_time_adapt_stab << " sec" << std::endl;
            std::cout << "****** Static condensation: " << newton_info.m_assembly_info.m_time_statcond << " sec" << std::endl;
            std::cout << "**** Solver time: " << newton_info.m_solve_info.m_time_solve << " sec" << std::endl;
            std::cout << "**** Postprocess time: " << newton_info.m_postprocess_info.m_time_post << " sec" << std::endl;
         }

         m_convergence = newton_solver.test_convergence();

         if(!m_convergence){
            if(step.level > m_rp.m_sublevel){
               std::cout << "***********************************************************" << std::endl;
               std::cout << "***** PROBLEM OF CONVERGENCE: We stop the calcul here *****" << std::endl;
               std::cout << "***********************************************************" << std::endl;
               break;
            }
            else
            {
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
            old_time = current_time;
            list_step.pop_front();
         }
      }

      if(m_convergence)
        newton_solver.save_solutions(m_solution_cells, m_solution_faces, m_solution_lagr, m_solution_data);

      ttot.toc();
      si.m_time_solver = ttot.to_double();
      return si;
   }

    bool test_convergence() const {return m_convergence;}

    template<typename AnalyticalSolution>
    scalar_type
    compute_l2_error(const AnalyticalSolution& as)
    {
        scalar_type err_dof = scalar_type{0.0};

        projector_type projk(m_bqd);

        size_t i = 0;

        for (auto& cl : m_msh)
        {
            auto x = m_solution_cells.at(i++);
            dynamic_vector<scalar_type> true_dof = projk.compute_cell(m_msh, cl, as);
            dynamic_vector<scalar_type> comp_dof = x.block(0,0,true_dof.size(), 1);
            dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
            err_dof += diff_dof.dot(projk.cell_mm * diff_dof);
        }

        return sqrt(err_dof);
    }

    size_t getDofs() {return total_dof_depl_static;}



    template<typename AnalyticalSolution>
    scalar_type
    compute_l2_gradient_error(const AnalyticalSolution& grad)
    {
       scalar_type err_dof = scalar_type{0.0};

       projector_type projk(m_bqd);

       gradrec_type gradrec(m_bqd);

       size_t i = 0;

       for (auto& cl : m_msh)
       {
          auto x = m_solution_data.at(i++);
          gradrec.compute(m_msh, cl, false);
          dynamic_vector<scalar_type> GTu = gradrec.oper()*x;
          dynamic_vector<scalar_type> true_dof = projk.compute_cell_grad(m_msh, cl, grad);
          dynamic_vector<scalar_type> comp_dof = GTu.block(0,0,true_dof.size(), 1);
          dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
          err_dof += diff_dof.dot(projk.grad_mm * diff_dof);
       }

       return sqrt(err_dof);
    }

    std::pair<scalar_type,scalar_type>
    compute_l2_error_annulus(const std::string& file_error)
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
          dynamic_vector<scalar_type> GTu = gradrec.oper()*x;

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
            auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, qp.point());
            auto GT_iqn = disk::compute_gradient_matrix_pt(GTu, gphi);
            auto FT_iqn = compute_FTensor(GT_iqn);
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

    std::array<scalar_type, 3>
    displacement_node(const size_t num_node)
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

    void
    compute_discontinuous_solution(const std::string& filename)
    {
       visu::Gmesh gmsh(DIM);
       auto storage = m_msh.backend_storage();

       std::vector<visu::Data> data; //create data (not used)
       std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point

       size_t cell_i(0);
       size_t nb_nodes(0);
       for (auto& cl : m_msh)
       {
          vector_dynamic x = m_solution_cells.at(cell_i++);
          auto cell_nodes = visu::cell_nodes(m_msh, cl);
          std::vector<visu::Node> new_nodes;
          for (size_t i = 0; i < cell_nodes.size(); i++)
          {
             nb_nodes++;
             auto point_ids = cell_nodes[i];
             auto pt = storage->points[point_ids];

             auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt);

             std::vector<scalar_type> depl(3, scalar_type{0});
             std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};

             visu::init_coor(pt, coor);
             visu::Node tmp_node(coor, nb_nodes, 0);
             new_nodes.push_back(tmp_node);
             gmsh.addNode(tmp_node);

             // plot magnitude at node
             for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_bqd.cell_degree()).size(); i += DIM)
                for(size_t j=0; j < DIM; j++)
                   depl[j] += phi.at(i+j)(j) * x(i+j); // a voir

               visu::Data datatmp(nb_nodes, depl);
               data.push_back(datatmp);
          }
          // add new element
          visu::add_element(gmsh, new_nodes);
       }

       visu::NodeData nodedata(3, 0.0, "depl_node", data, subdata); // create and init a nodedata view

       nodedata.saveNodeData(filename, gmsh); // save the view
    }


    void
    compute_conforme_solution(const std::string& filename)
    {
       visu::Gmesh gmsh = visu::convertMesh(m_msh);
       auto storage = m_msh.backend_storage();
       std::vector<double> vzero(3,0.0);
       size_t nb_nodes(gmsh.getNumberofNodes());

       //first(number of data at this node), second(cumulated value)
       std::vector<std::pair<size_t, std::vector<scalar_type> > > value(nb_nodes, std::make_pair(0, vzero));

       size_t cell_i(0);
       for (auto& cl : m_msh)
       {
          vector_dynamic x = m_solution_cells.at(cell_i++);
          auto cell_nodes = visu::cell_nodes(m_msh, cl);
          std::vector<visu::Node> new_nodes;
          for (size_t i = 0; i < cell_nodes.size(); i++)
          {
             nb_nodes++;
             auto point_ids = cell_nodes[i];
             auto pt = storage->points[point_ids];

             auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt);

             std::vector<scalar_type> depl(3, scalar_type{0});

             // plot magnitude at node
             for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_bqd.cell_degree()).size(); i += DIM)
                for(size_t j=0; j < DIM; j++)
                   depl[j] += phi.at(i+j)(j) * x(i+j); // a voir

             value[point_ids].first +=1;
             for(size_t j=0; j < DIM; j++)
               value[point_ids].second[j] += depl[j];
          }
       }

       std::vector<visu::Data> data; //create data
       std::vector<visu::SubData> subdata; //create subdata
       data.reserve(nb_nodes); // data has a size of nb_node

      for(size_t  i_node = 0; i_node < value.size(); i_node++){
         std::vector<double> tmp_value(3,0.0);
         for(size_t j=0; j < DIM; j++)
            tmp_value[j] = value[i_node].second[j]/ double(value[i_node].first);

         visu::Data tmp_data(i_node + 1, tmp_value);
         data.push_back(tmp_data); //add data
      }

       visu::NodeData nodedata(3, 0.0, "depl_node", data, subdata); // create and init a nodedata view

       nodedata.saveNodeData(filename, gmsh); // save the view
    }

    void
    compute_deformed(const std::string& filename)
    {
       visu::Gmesh gmsh(DIM);
       auto storage = m_msh.backend_storage();

       size_t cell_i(0);
       size_t nb_nodes(0);
       for (auto& cl : m_msh)
       {
          vector_dynamic x = m_solution_cells.at(cell_i++);
          auto cell_nodes = visu::cell_nodes(m_msh, cl);
          std::vector<visu::Node> new_nodes;
          for (size_t i = 0; i < cell_nodes.size(); i++)
          {
             nb_nodes++;
             auto point_ids = cell_nodes[i];
             auto pt = storage->points[point_ids];

             auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt);

             std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};
             std::array<double, 3> depl = {double{0.0}, double{0.0}, double{0.0}};

             visu::init_coor(pt, coor);
             for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_bqd.cell_degree()).size(); i += DIM)
                for(size_t j=0; j < DIM; j++)
                   depl[j] += phi.at(i+j)(j) * x(i+j); // a voir

            for(size_t j=0; j < DIM; j++)
               coor[j] += depl[j];

            visu::Node tmp_node(coor, nb_nodes, 0);
            new_nodes.push_back(tmp_node);
            gmsh.addNode(tmp_node);

          }
          // add new element
          visu::add_element(gmsh, new_nodes);
       }
       gmsh.writeGmesh(filename, 2);

    }

    void
    compute_discontinuous_PK1(const std::string& filename)
    {
       visu::Gmesh gmsh(DIM);
       auto storage = m_msh.backend_storage();

       gradrec_type gradrec(m_bqd);
       NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);


       std::vector<visu::Data> data; //create data (not used)
       std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point

       size_t cell_i(0);
       size_t nb_nodes(0);
       for (auto& cl : m_msh)
       {
          gradrec.compute(m_msh, cl, false);
          const vector_dynamic GT_uTF = gradrec.oper() * m_solution_data.at(cell_i);

          auto cell_nodes = visu::cell_nodes(m_msh, cl);
          std::vector<visu::Node> new_nodes;
          for (size_t i = 0; i < cell_nodes.size(); i++)
          {
             nb_nodes++;
             auto point_ids = cell_nodes[i];
             auto pt = storage->points[point_ids];

             auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, pt);

             auto GT_iqn = disk::compute_gradient_matrix_pt(GT_uTF, gphi);
             auto FT_iqn = compute_FTensor(GT_iqn);

             auto PK1= law.compute_PK1(FT_iqn);

             std::vector<scalar_type> PK1_pt(1, PK1.norm());
             std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};

             visu::init_coor(pt, coor);
             visu::Node tmp_node(coor, nb_nodes, 0);
             new_nodes.push_back(tmp_node);
             gmsh.addNode(tmp_node);


               visu::Data datatmp(nb_nodes, PK1_pt);
               data.push_back(datatmp);
          }
          // add new element
          visu::add_element(gmsh, new_nodes);
          cell_i++;
       }

       visu::NodeData nodedata(1, 0.0, "PK1", data, subdata); // create and init a nodedata view

       nodedata.saveNodeData(filename, gmsh); // save the view
    }

    //compute PK in cylindrical base
    void
    compute_discontinuous_Prr(const std::string& filename, const std::string compo = "Prr")
    {
       visu::Gmesh gmsh(DIM);
       auto storage = m_msh.backend_storage();

       gradrec_type gradrec(m_bqd);
       NeoHookeanLaw<scalar_type>  law(m_elas_param.mu, m_elas_param.lambda, m_elas_param.type_law);


       std::vector<visu::Data> data; //create data (not used)
       std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point

       size_t cell_i(0);
       size_t nb_nodes(0);
       for (auto& cl : m_msh)
       {
          gradrec.compute(m_msh, cl, false);
          const vector_dynamic GT_uTF = gradrec.oper() * m_solution_data.at(cell_i);

          auto cell_nodes = visu::cell_nodes(m_msh, cl);
          std::vector<visu::Node> new_nodes;
          for (size_t i = 0; i < cell_nodes.size(); i++)
          {
             nb_nodes++;
             auto point_ids = cell_nodes[i];
             auto pt = storage->points[point_ids];

             auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, pt);

             auto GT_iqn = disk::compute_gradient_matrix_pt(GT_uTF, gphi);
             auto FT_iqn = compute_FTensor(GT_iqn);

             auto PK1= law.compute_PK1(FT_iqn);

             vector_dynamic er;
             er.resize(2);
             er(0) = pt.x(); er(1) = pt.y();
             er /= er.norm();
             vector_dynamic eO;
             eO.resize(2);
             eO(0) = -er(1); eO(1) = er(0);

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


    template<typename AnalyticalSolution>
    void
    plot_l2error_at_gausspoint(const std::string& filename, const AnalyticalSolution& as)
    {
       visu::Gmesh msh; //creta a mesh

       std::vector<visu::Data> data; //create data (not used)
       std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point
       size_t nb_node =  msh.getNumberofNodes();

       size_t cell_i = 0;
       for (auto& cl : m_msh)
       {
          vector_dynamic x = m_solution_cells.at(cell_i++);
          auto qps = m_bqd.cell_quadrature.integrate(m_msh, cl);
          for (auto& qp : qps)
          {

             auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, qp.point());

             std::array<double, 3> depl = {double{0.0}, double{0.0}, double{0.0}};

             for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_bqd.m_cell_degree()).size(); i += DIM)
                for(size_t j = 0; j < DIM; j++)
                   depl[j] += phi.at(i+j)(j) * x(i+j); // a voir

             auto true_depl = as(qp.point()); // a voir et projeté

             nb_node += 1;
             visu::Node snode = visu::convertPoint(qp.point(), nb_node); //create a node at gauss point
             std::vector<double> value(1, 0.0); // save the solution at gauss point
             for(size_t i = 0; i < DIM; i++)
                value[0] += std::pow(depl[i] - true_depl(i), 2.0);

             value[0] = sqrt(value[0]);

             visu::SubData sdata(value, snode);
             subdata.push_back(sdata); // add subdata


          }
       }

       visu::NodeData nodedata(1, 0.0, "error_depl", data, subdata); // create and init a nodedata view

       nodedata.saveNodeData(filename, msh); // save the view
    }

    void
    plot_displacement_at_gausspoint(const std::string& filename)
    {
       visu::Gmesh msh; //creta a mesh

       std::vector<visu::Data> data; //create data (not used)
       std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point
       size_t nb_node =  msh.getNumberofNodes();

       size_t cell_i = 0;
       for (auto& cl : m_msh)
       {
          vector_dynamic x = m_solution_cells.at(cell_i++);
          auto qps = m_bqd.cell_quadrature.integrate(m_msh, cl);
          for (auto& qp : qps)
          {

             auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, qp.point());

             std::vector<double> depl(3, double{0.0});

             for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_bqd.cell_degree()).size(); i += DIM)
                for(size_t j = 0; j < DIM; j++)
                   depl[j] += phi.at(i+j)(j) * x(i+j); // a voir


             nb_node += 1;
             visu::Node snode = visu::convertPoint(qp.point(), nb_node); //create a node at gauss point
             visu::SubData sdata(depl, snode);
             subdata.push_back(sdata); // add subdata


          }
       }

       visu::NodeData nodedata(3, 0.0, "depl", data, subdata); // create and init a nodedata view

       nodedata.saveNodeData(filename, msh); // save the view
    }


    void
    plot_J_at_gausspoint(const std::string& filename)
    {
       visu::Gmesh msh; //creta a mesh

       std::vector<visu::Data> data; //create data (not used)
       std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point
       size_t nb_node =  msh.getNumberofNodes();

       gradrec_type gradrec(m_bqd);

       size_t cell_i = 0;
       for (auto& cl : m_msh)
       {
          vector_dynamic x = m_solution_data.at(cell_i++);
          gradrec.compute(m_msh, cl, false);
          vector_dynamic GTu = gradrec.oper()*x;

          auto qps = m_bqd.grad_quadrature.integrate(m_msh, cl);
          for (auto& qp : qps)
          {

             auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, qp.point());

             // Compute local gradient and norm
             auto GT_iqn = disk::compute_gradient_matrix_pt(GTu, gphi);
             auto FT_iqn = compute_FTensor(GT_iqn);

             std::vector<double> J(1, FT_iqn.determinant());

             nb_node += 1;
             visu::Node snode = visu::convertPoint(qp.point(), nb_node); //create a node at gauss point
             visu::SubData sdata(J, snode);
             subdata.push_back(sdata); // add subdata


          }
       }

       visu::NodeData nodedata(1, 0.0, "Jacobian", data, subdata); // create and init a nodedata view

       nodedata.saveNodeData(filename, msh); // save the view
    }


    void
    plot_J(const std::string& filename)
    {
       visu::Gmesh gmsh(DIM);
       auto storage = m_msh.backend_storage();

       std::vector<visu::Data> data; //create data (not used)
       std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point

       gradrec_type gradrec(m_bqd);

       size_t cell_i(0);
       size_t nb_nodes(0);
       for (auto& cl : m_msh)
       {
          vector_dynamic x = m_solution_data.at(cell_i++);
          gradrec.compute(m_msh, cl, false);
          vector_dynamic GTu = gradrec.oper()*x;
          auto cell_nodes = visu::cell_nodes(m_msh, cl);
          std::vector<visu::Node> new_nodes;
          for (size_t i = 0; i < cell_nodes.size(); i++)
          {
             nb_nodes++;
             auto point_ids = cell_nodes[i];
             auto pt = storage->points[point_ids];

             std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};

             visu::init_coor(pt, coor);
             visu::Node tmp_node(coor, nb_nodes, 0);
             new_nodes.push_back(tmp_node);
             gmsh.addNode(tmp_node);


             auto gphi = m_bqd.grad_basis.eval_functions(m_msh, cl, pt);

             // Compute local gradient and norm
             auto GT_iqn = disk::compute_gradient_matrix_pt(GTu, gphi);
             auto FT_iqn = compute_FTensor(GT_iqn);

             std::vector<double> J(1, FT_iqn.determinant());

             visu::Data datatmp(nb_nodes, J);
             data.push_back(datatmp);

          }
          // add new element
          visu::add_element(gmsh, new_nodes);
       }

       visu::NodeData nodedata(1, 0.0, "J_node", data, subdata); // create and init a nodedata view

       nodedata.saveNodeData(filename, gmsh); // save the view
    }

};
