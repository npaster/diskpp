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
#include "hho_nl.hpp"
#include "bases/bases_utils.hpp"
#include "NewtonSolver/newton_solver.hpp"


#include "../exemple_visualisation/visualisation/gmshDisk.hpp"
#include "../exemple_visualisation/visualisation/gmshConvertMesh.hpp"
#include "../exemple_visualisation/visualisation/gmshDeformed.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>


struct offline_info
{
    double  time_offline;
};

struct solve_info
{
    double  time_solver;
};


template<template<typename, size_t , typename> class Mesh,
typename T, size_t DIM, typename Storage, typename Point>
class NL_elasticity_solver
{
   typedef Mesh<T, DIM, Storage>                      mesh_type;
   typedef typename mesh_type::scalar_type            scalar_type;
   typedef typename mesh_type::cell                   cell_type;
   typedef typename mesh_type::face                   face_type;

   typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
   typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;

   typedef disk::scaled_monomial_vector_basis<mesh_type, cell_type>    cell_basis_type;
   typedef disk::scaled_monomial_vector_basis<mesh_type, face_type>    face_basis_type;

   typedef dynamic_matrix<scalar_type>         matrix_dynamic;
   typedef dynamic_vector<scalar_type>         vector_dynamic;

   size_t m_cell_degree, m_face_degree, m_degree;

   const mesh_type& m_msh;

   std::vector<vector_dynamic>                    m_solution_data;
   std::vector<matrix_dynamic>                    m_data_offline;
   std::vector<vector_dynamic>         m_solution_cells, m_solution_faces, m_solution_lagr;

   bool m_verbose, m_convergence;

public:
   NL_elasticity_solver(const mesh_type& msh, size_t degree, int l = 0)
   : m_msh(msh), m_verbose(false), m_convergence(false)
   {
      if ( l < -1 or l > 1)
      {
         std::cout << "'l' should be -1, 0 or 1. Reverting to 0." << std::endl;
         l = 0;
      }

      if (degree == 0 && l == -1)
      {
         std::cout << "'l' should be 0 or 1. Reverting to 0." << std::endl;
         l = 0;
      }

      m_cell_degree = degree + l;
      m_face_degree = degree;
      m_degree = degree;

   }

   bool    verbose(void) const     { return m_verbose; }
   void    verbose(bool v)         { m_verbose = v; }

   offline_info
   compute_offline()
   {
      disk::gradient_reconstruction_elas< mesh_type,
                                          cell_basis_type,
                                          cell_quadrature_type,
                                          face_basis_type,
                                          face_quadrature_type> gradrec(m_degree);


//       disk::diffusion_like_stabilization< mesh_type,
//                                           cell_basis_type,
//                                           cell_quadrature_type,
//                                           face_basis_type,
//                                           face_quadrature_type> stab(m_degree);

      m_data_offline.reserve(m_msh.cells_size());

      offline_info ai;
      bzero(&ai, sizeof(ai));

      timecounter tc;

      tc.tic();
      for (auto& cl : m_msh)
      {
         gradrec.compute(m_msh, cl);
//          stab.compute(m_msh, cl, gradrec.oper);
         //dynamic_matrix<scalar_type> loc = gradrec.oper;// + stab.data;
         m_data_offline.push_back(gradrec.oper);
      }
      tc.toc();
      ai.time_offline += tc.to_double();

      assert(m_data_offline.size() == m_msh.cells_size());

      return ai;
   }

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

      cell_basis_type cell_basis = cell_basis_type(m_degree);
      face_basis_type face_basis = face_basis_type(m_degree);

      const size_t num_cell_dofs = cell_basis.size();
      const size_t num_face_dofs = face_basis.size();
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
         std::cout << "** Numbers of dofs: " << total_dof + total_lagr  << std::endl;
         std::cout << "** including " << total_lagr << " Lagrange multipliers"  << std::endl;
         std::cout << "** After static condensation: "  << std::endl;
         std::cout << "** Numbers of dofs: " << m_msh.faces_size() * num_face_dofs + total_lagr  << std::endl;
         std::cout << "** including " << total_lagr << " Lagrange multipliers"  << std::endl;
      }
      //provisoire

//       for(size_t i = 0; i < m_msh.cells_size(); i++)
//       {
//          m_solution_data[i].setConstant(0.001);
//          m_solution_cells[i].setConstant(0.001);
//       }
//
//       for(size_t i = 0; i < m_msh.faces_size(); i++)
//       {
//          m_solution_faces[i].setConstant(0.001);
//       }
//
//       for(size_t i = 0; i < m_msh.boundary_faces_size(); i++){
//          m_solution_lagr[i].setConstant(0.001);
//      }
   }

   template<typename LoadFunction, typename BoundaryConditionFunction>
   solve_info
   compute(const LoadFunction& lf, const BoundaryConditionFunction& bcf,
      size_t n_time_step = 1)
   {

      solve_info ai;
      bzero(&ai, sizeof(ai));

      timecounter tc;

      NewtonRaphson_solver_elasticity<Mesh<T,DIM,Storage>> newton_solver(m_msh, m_degree, 0);

      newton_solver.initialize(m_solution_cells, m_solution_faces,
                                 m_solution_lagr, m_solution_data);

      const scalar_type delta_t = 1.0/n_time_step;

      scalar_type time = 0.0;

      for (size_t n = 0; n < n_time_step; n++)
      {
         time += delta_t;
         tc.tic();
         if(m_verbose){
            std::cout << "----------------------------------------------" << std::endl;
            std::cout << "************** Time step " << n+1 << "/" << n_time_step << " *****************" << std::endl;
         }

         auto rlf = [&lf, &time](const Point& p) -> auto {
             return disk::mm_prod(time,lf(p));

         };

         auto rbcf = [&bcf, &time](const Point& p) -> auto {
             return disk::mm_prod(time,bcf(p));
         };

         auto newton_info = newton_solver.compute(rlf, rbcf, m_data_offline);

         tc.toc();
         ai.time_solver += tc.to_double();

         if(m_verbose){
            std::cout << "** Time in this step " << tc.to_double() << " sec" << std::endl;
            std::cout << "**** Assembly time: " << newton_info.time_assembly << " sec" << std::endl;
            std::cout << "****** Gradient reconstruction: " << newton_info.time_gradrec << " sec" << std::endl;
            std::cout << "****** Stabilisation: " << newton_info.time_stab << " sec" << std::endl;
            std::cout << "****** Elementary computation: " << newton_info.time_elem << " sec" << std::endl;
            std::cout << "****** Static condensation: " << newton_info.time_statcond << " sec" << std::endl;
            std::cout << "**** Solver time: " << newton_info.time_solve << " sec" << std::endl;
            std::cout << "**** Postprocess time: " << newton_info.time_post << " sec" << std::endl;
         }

         m_convergence = newton_solver.test_convergence();

         if(!m_convergence){
             std::cout << "***********************************************************" << std::endl;
             std::cout << "***** PROBLEM OF CONVERGENCE: We stop the calcul here *****" << std::endl;
             std::cout << "***********************************************************" << std::endl;
             break;
         }
      }

      if(m_convergence)
        newton_solver.save_solutions(m_solution_cells, m_solution_faces,
                                 m_solution_lagr, m_solution_data);

      return ai;
   }

    bool test_convergence() const {return m_convergence;}

    template<typename AnalyticalSolution>
    scalar_type
    compute_l2_error(const AnalyticalSolution& as)
    {
        scalar_type err_dof = scalar_type{0.0};

        disk::projector_elas<mesh_type, cell_basis_type, cell_quadrature_type,
                        face_basis_type, face_quadrature_type> projk(m_degree);

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


void
    plot_solution_at_gausspoint(const std::string& filename)
    {
       std::cout << "Compute solution at Gauss points" << std::endl;
       visu::Gmesh msh; //creta a mesh

        std::vector<visu::Data> data; //create data (not used)
        std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point
        size_t nb_node =  msh.getNumberofNodes();

        size_t dim = m_msh.dimension;

        cell_basis_type cell_basis          = cell_basis_type(m_degree);
        cell_quadrature_type cell_quadrature     = cell_quadrature_type(m_degree);

        size_t cell_i = 0;
        for (auto& cl : m_msh)
        {
            vector_dynamic x = m_solution_cells.at(cell_i++);
            auto qps = cell_quadrature.integrate(m_msh, cl);
            for (auto& qp : qps)
            {

                auto phi = cell_basis.eval_functions(m_msh, cl, qp.point());

                vector_dynamic pot = vector_dynamic::Zero(dim);
                for (size_t i = 0; i < cell_basis.range(0, m_cell_degree).size(); i+=dim)
                    for(size_t j=0; j<dim; j++)
                            pot(j) += phi.at(i+j)(j) * x(i+j);


               nb_node += 1;
               visu::Node snode = visu::convertPoint(qp.point(), nb_node); //create a node at gauss point
               std::vector<double> value(3, 0.0);
               for(size_t j=0; j<dim; j++)
                   value.at(j) =  {double(pot(j))}; // save the solution at gauss point
               visu::SubData sdata(value, snode);
               subdata.push_back(sdata); // add subdata


            }
        }

        visu::NodeData nodedata(3, 0.0, "sol_gp", data, subdata); // create and init a nodedata view

        nodedata.saveNodeData(filename, msh); // save the view
    }


   template<typename AnalyticalSolution>
   void
   plot_l2error_at_gausspoint(const std::string& filename, const AnalyticalSolution& as)
   {
      std::cout << "Compute L2 error at Gauss points" << std::endl;
      visu::Gmesh msh; //creta a mesh

      std::vector<visu::Data> data; //create data (not used)
      std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point
      size_t nb_node =  msh.getNumberofNodes();

      size_t dim = m_msh.dimension;

      cell_basis_type cell_basis          = cell_basis_type(m_degree);
      cell_quadrature_type cell_quadrature     = cell_quadrature_type(m_degree);

      size_t cell_i = 0;
      for (auto& cl : m_msh)
      {
            vector_dynamic x = m_solution_cells.at(cell_i++);
            auto qps = cell_quadrature.integrate(m_msh, cl);
            for (auto& qp : qps)
            {

               auto phi = cell_basis.eval_functions(m_msh, cl, qp.point());

               vector_dynamic pot = vector_dynamic::Zero(dim);
               for (size_t i = 0; i < cell_basis.range(0, m_cell_degree).size(); i+=dim)
                  for(size_t j=0; j<dim; j++)
                           pot(j) += phi.at(i+j)(j) * x(i+j);

               auto true_pot = as(qp.point()); // a voir et projeté

               nb_node += 1;
               visu::Node snode = visu::convertPoint(qp.point(), nb_node); //create a node at gauss point
               std::vector<double> value(1,0.0);
               double l2error = 0.0;
               for(size_t j=0; j<dim; j++)
                  l2error +=  double((pot(j)-true_pot(j))*(pot(j)-true_pot(j))); // save the solution at gauss point

               value[0] = sqrt(l2error);
               visu::SubData sdata(value, snode);
               subdata.push_back(sdata); // add subdata


            }
      }

      visu::NodeData nodedata(1, 0.0, "l2error_gp", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, msh); // save the view
   }

   void
    saveMesh2D(const std::string& filename)
    {
       visu::Gmesh gmsh = visu::convertMesh2D(m_msh);
       gmsh.writeGmesh(filename, 2);
    }

     void
     saveMesh3D(const std::string& filename)
     {
        visu::Gmesh gmsh = visu::convertMesh3D(m_msh);
        gmsh.writeGmesh(filename, 2);
     }

   void
   compute_deformed(const std::string& filename)
   {
      visu::Gmesh gmsh(DIM);
      auto storage = m_msh.backend_storage();

      cell_basis_type cell_basis          = cell_basis_type(m_degree);
      cell_quadrature_type cell_quadrature     = cell_quadrature_type(m_degree);

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

            auto phi = cell_basis.eval_functions(m_msh, cl, pt);

            std::vector<scalar_type> coor(3, scalar_type{0});
            std::vector<scalar_type> depl(3, scalar_type{0});

            visu::init_coor(pt, coor);
            for (size_t i = 0; i < cell_basis.range(0, m_cell_degree).size(); i += DIM)
               for(size_t j=0; j < DIM; j++)
                  depl[j] += phi.at(i+j)(j) * x(i+j); // a voir

            for(size_t j=0; j < DIM; j++)
               coor[j] += depl[j];

            visu::Node tmp_node(coor, nb_nodes, 0);
            new_nodes.push_back(tmp_node);
            gmsh.addNode(tmp_node);

         }
         // add new element
         visu::add_element<DIM>(gmsh, new_nodes);
      }
      gmsh.writeGmesh(filename, DIM);

   }


   void
   compute_discontinuous_solution(const std::string& filename)
   {
      visu::Gmesh gmsh(DIM);
      auto storage = m_msh.backend_storage();

      cell_basis_type cell_basis          = cell_basis_type(m_degree);
      cell_quadrature_type cell_quadrature     = cell_quadrature_type(m_degree);

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

            auto phi = cell_basis.eval_functions(m_msh, cl, pt);

            std::vector<scalar_type> coor(3, scalar_type{0});
            std::vector<scalar_type> depl(3, scalar_type{0});

            visu::init_coor(pt, coor);
            visu::Node tmp_node(coor, nb_nodes, 0);
            new_nodes.push_back(tmp_node);
            gmsh.addNode(tmp_node);

            // plot magnitude at node
            for (size_t i = 0; i < cell_basis.range(0, m_cell_degree).size(); i += DIM)
               for(size_t j=0; j < DIM; j++)
                  depl[j] += phi.at(i+j)(j) * x(i+j); // a voir

            visu::Data datatmp(nb_nodes, depl);
            data.push_back(datatmp);
         }
         // add new element
         visu::add_element<DIM>(gmsh, new_nodes);
      }

      visu::NodeData nodedata(3, 0.0, "depl_node", data, subdata); // create and init a nodedata view

      nodedata.saveNodeData(filename, gmsh); // save the view
   }

//     void
//     compute_deformed1D(const std::string& filename)
//     {
//        visu::Gmesh gmsh(1);
//        auto storage = m_msh.backend_storage();
//
//        cell_basis_type cell_basis          = cell_basis_type(m_degree);
//        cell_quadrature_type cell_quadrature     = cell_quadrature_type(m_degree);
//
//        size_t cell_i(0);
//        size_t nb_nodes(0);
//        for (auto& cl : m_msh)
//        {
//           vector_dynamic x = m_solution_cells.at(cell_i++);
//           auto cell_nodes = visu::cell_nodes(m_msh, cl);
//           std::vector<visu::Node> new_nodes;
//           for (size_t i = 0; i < cell_nodes.size(); i++)
//           {
//              nb_nodes++;
//              auto point_ids = cell_nodes[i];
//              auto pt = storage->points[point_ids];
//
//              auto phi = cell_basis.eval_functions(m_msh, cl, pt);
//
//              std::vector<scalar_type> coor(3, scalar_type{0});
//              for (size_t i = 0; i < cell_basis.range(0, m_cell_degree).size(); i++)
//                 coor[0] += phi.at(i) * x(i); // a voir
//
//                 visu::Node tmp_node(coor, nb_nodes, 0);
//              new_nodes.push_back(tmp_node);
//              gmsh.addNode(tmp_node);
//           }
//           // add new element
//        }
//        gmsh.writeGmesh(filename, 2);
//     }
};
