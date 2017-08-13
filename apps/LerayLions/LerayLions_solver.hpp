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

#include <sstream>

#include "../../config.h"

#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "hho/hho.hpp"
#include "hho/hho_vector.hpp"
#include "bases/bases_utils.hpp"
#include "NewtonSolver/newton_solver.hpp"
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
class leraylions_solver
{
   typedef Mesh<T, DIM, Storage>                      mesh_type;
   typedef typename mesh_type::scalar_type            scalar_type;
   typedef ParamRun<T>                                param_type;

   typedef disk::basis_quadrature_data_full<mesh_type, disk::scaled_monomial_scalar_basis,
                                                      disk::scaled_monomial_vector_basis,
                                                      disk::quadrature> bqdata_type;

   typedef dynamic_matrix<scalar_type>         matrix_dynamic;
   typedef dynamic_vector<scalar_type>         vector_dynamic;

   typedef disk::gradient_reconstruction_full_bq<bqdata_type>     gradrec_type;
   typedef disk::projector_bq<bqdata_type>                        projector_type;

   bqdata_type     m_bqd;

   const mesh_type& m_msh;

   std::vector<vector_dynamic>         m_solution_data;
   std::vector<vector_dynamic>         m_solution_cells, m_solution_faces, m_solution_lagr;

   bool m_verbose, m_convergence;

   scalar_type m_leray_param;
   param_type m_rp;

   size_t total_dof_depl_static;

public:
   leraylions_solver(const mesh_type& msh, const param_type& rp, const scalar_type leray_param)
   : m_msh(msh), m_verbose(rp.m_verbose), m_convergence(false), m_leray_param(leray_param), m_rp(rp)
   {
      size_t face_degree = rp.m_face_degree;
      if(face_degree <= 0)
      {
         std::cout << "'face_degree' should be > 0. Reverting to 1." << std::endl;
         face_degree = 1;
      }

      m_rp.m_face_degree = face_degree;

      size_t cell_degree = rp.m_cell_degree;
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

      m_solution_data.reserve(m_msh.cells_size());
      m_solution_cells.reserve(m_msh.cells_size());
      m_solution_faces.reserve(m_msh.faces_size());
      m_solution_lagr.reserve(m_msh.boundary_faces_size());


      const size_t num_cell_dofs = (m_bqd.cell_basis.range(0, m_bqd.cell_degree())).size();
      const size_t num_face_dofs = m_bqd.face_basis.size();
      const size_t total_dof = m_msh.cells_size() * num_cell_dofs + m_msh.faces_size() * num_face_dofs;
      const size_t total_lagr = m_msh.boundary_faces_size() * num_face_dofs;

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

      for(size_t i = 0; i < m_msh.boundary_faces_size(); i++){
         m_solution_lagr.push_back(vector_dynamic::Zero(num_face_dofs));
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

      for(size_t i = 0; i < m_msh.cells_size(); i++)
      {
         m_solution_data[i].setConstant(0.01);
         m_solution_cells[i].setConstant(0.01);
      }

      for(size_t i = 0; i < m_msh.faces_size(); i++)
      {
         m_solution_faces[i].setConstant(0.01);
      }

      for(size_t i = 0; i < m_msh.boundary_faces_size(); i++){
         m_solution_lagr[i].setConstant(0.01);
     }
   }

   template<typename LoadFunction, typename BoundaryConditionFunction>
   SolverInfo
   compute(const LoadFunction& lf, const BoundaryConditionFunction& bcf)
   {

      SolverInfo si;
      timecounter ttot;
      ttot.tic();

      //time step
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


      //Newton solver
      NewtonRaphson_solver_leraylions<bqdata_type> newton_solver(m_msh, m_bqd, m_leray_param);

      newton_solver.initialize(m_solution_cells, m_solution_faces, m_solution_lagr, m_solution_data);

      newton_solver.verbose(m_verbose);


      //loading
      while(!list_step.empty())
      {
         current_step += 1;
         time_step step = list_step.front();
         const scalar_type current_time = step.time;
         if(m_verbose){
            std::cout << "---------------------------------------------------------------------------------------------------------------------" << std::endl;
            std::cout << "****************************** Time : " << current_time << " sec (step: " << current_step
            << "/" << total_step << ", sublevel: " << step.level << " ) *******************************|" << std::endl;
         }

         auto rlf = [&lf, &current_time](const Point& p) -> auto {
            return disk::mm_prod(current_time, lf(p));

         };

         auto rbcf = [&bcf, &current_time](const Point& p) -> auto {
            return disk::mm_prod(current_time, bcf(p));
         };

         //correction
         NewtonSolverInfo newton_info = newton_solver.compute(rlf, rbcf, m_rp.m_epsilon, m_rp.m_iter_max);
         si.updateInfo(newton_info);

         if(m_verbose){
            std::cout << "** Time in this step " << newton_info.m_time_newton << " sec" << std::endl;
            std::cout << "**** Assembly time: " << newton_info.m_assembly_info.m_time_assembly << " sec" << std::endl;
            std::cout << "****** Gradient reconstruction: " << newton_info.m_assembly_info.m_time_gradrec << " sec" << std::endl;
            std::cout << "****** Elementary computation: " << newton_info.m_assembly_info.m_time_elem << " sec" << std::endl;
            std::cout << "       *** Behavior computation: " << newton_info.m_assembly_info.m_time_law << " sec" << std::endl;
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
          gradrec.compute_optim(m_msh, cl,false);
          dynamic_vector<scalar_type> GTu = gradrec.oper()*x;
          dynamic_vector<scalar_type> true_dof = projk.compute_cell_grad(m_msh, cl, grad);
          dynamic_vector<scalar_type> comp_dof = GTu.block(0,0,true_dof.size(), 1);
          dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
          err_dof += diff_dof.dot(projk.grad_mm * diff_dof);
       }

       return sqrt(err_dof);
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

             std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};

             visu::init_coor(pt, coor);
             visu::Node tmp_node(coor, nb_nodes, 0);
             new_nodes.push_back(tmp_node);
             gmsh.addNode(tmp_node);

             // plot magnitude at node
             std::vector<double> sol(1, x.dot(phi));

             visu::Data datatmp(nb_nodes, sol);
             data.push_back(datatmp);
          }
          // add new element
          visu::add_element(gmsh, new_nodes);
       }

       visu::NodeData nodedata(1, 0.0, "sol", data, subdata); // create and init a nodedata view

       nodedata.saveNodeData(filename, gmsh); // save the view
    }


    void
    compute_conforme_solution(const std::string& filename)
    {
       visu::Gmesh gmsh = visu::convertMesh(m_msh);
       auto storage = m_msh.backend_storage();
       size_t nb_nodes(gmsh.getNumberofNodes());

       //first(number of data at this node), second(cumulated value)
       std::vector<std::pair<size_t, scalar_type>> value(nb_nodes, std::make_pair(0, 0.0));

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

             // plot magnitude at node
             scalar_type sol = phi.dot(x);

             value[point_ids].first +=1;
             value[point_ids].second += sol;
          }
       }

       std::vector<visu::Data> data; //create data
       std::vector<visu::SubData> subdata; //create subdata
       data.reserve(nb_nodes); // data has a size of nb_node

      for(size_t  i_node = 0; i_node < value.size(); i_node++){
         std::vector<double> tmp_value(1, value[i_node].second/ double(value[i_node].first));

         visu::Data tmp_data(i_node + 1, tmp_value);
         data.push_back(tmp_data); //add data
      }

       visu::NodeData nodedata(1, 0.0, "sol", data, subdata); // create and init a nodedata view

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
             scalar_type sol = phi.dot(x);

             auto true_depl = as(qp.point()); // a voir et projeté

             nb_node += 1;
             visu::Node snode = visu::convertPoint(qp.point(), nb_node); //create a node at gauss point
             std::vector<double> value(1, std::abs(sol - true_depl)); // save the solution at gauss point
             visu::SubData sdata(value, snode);
             subdata.push_back(sdata); // add subdata
          }
       }

       visu::NodeData nodedata(1, 0.0, "error_sol", data, subdata); // create and init a nodedata view

       nodedata.saveNodeData(filename, msh); // save the view
    }

    void
    plot_solution_at_gausspoint(const std::string& filename)
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
             std::vector<double> sol(1, phi.dot(x));

             nb_node += 1;
             visu::Node snode = visu::convertPoint(qp.point(), nb_node); //create a node at gauss point
             visu::SubData sdata(sol, snode);
             subdata.push_back(sdata); // add subdata
          }
       }

       visu::NodeData nodedata(1, 0.0, "sol", data, subdata); // create and init a nodedata view

       nodedata.saveNodeData(filename, msh); // save the view
    }


    void
    compute_deformed(const std::string& filename)
    {
       if(DIM >= 3)
          std::cout << "Compute deformed only in 1D or 2D" << '\n';
       else {
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

                std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};

                visu::init_coor(pt, coor);
                visu::Node tmp_node(coor, nb_nodes, 0);
                new_nodes.push_back(tmp_node);
                gmsh.addNode(tmp_node);

                // plot magnitude at node
                scalar_type pot = x.dot(phi);

                std::vector<double> value(3, 0.0);
                value[DIM+1 -1] = {double(pot)};
                visu::Data datatmp(nb_nodes, value);

                data.push_back(datatmp);
             }
             // add new element
             visu::add_element(gmsh, new_nodes);
          }

          visu::NodeData nodedata(3, 0.0, "sol_node_deformed", data, subdata); // create and init a nodedata view

          nodedata.saveNodeData(filename, gmsh); // save the view
       }
    }

};
