/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
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
#include <sstream>

#include <list>
#include <vector>

#include "bases/bases.hpp"
#include "methods/hho"
#include "quadratures/quadratures.hpp"

#include "solvers/solver.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

  template<typename Mesh>
  class diffusion_solver
{
    typedef Mesh                                       mesh_type;
    typedef typename mesh_type::coordinate_type        scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;


    typename disk::hho_degree_info m_hdi;

    const mesh_type& m_msh;

    std::vector<vector_type> m_postprocess_data;

    bool m_verbose;


public:
   diffusion_solver(const mesh_type& msh, size_t degree, int l = 0)
   : m_msh(msh), m_verbose(false)
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

      m_hdi         = disk::hho_degree_info(degree + l, degree);
   }

   bool    verbose(void) const     { return m_verbose; }
   void    verbose(bool v)         { m_verbose = v; }

   template<typename LoadFunction, typename BoundaryConditionFunction>
   void
   solve(const LoadFunction& lf, const BoundaryConditionFunction& bcf)
   {
       auto assembler = make_diffusion_assembler(m_msh, m_hdi);

       for (auto& cl : m_msh)
       {
           const auto gr = make_scalar_hho_laplacian(m_msh, cl, m_hdi);

           const auto stab = make_scalar_hho_stabilization(m_msh, cl, gr.first, m_hdi);

           auto cb  = make_scalar_monomial_basis(m_msh, cl, m_hdi.cell_degree());
           auto rhs = make_rhs(m_msh, cl, cb, lf);
           matrix_type A   = gr.second + stab;
           auto sc  = make_scalar_static_condensation(m_msh, cl, m_hdi, A, rhs);

           assembler.assemble(m_msh, cl, sc.first, sc.second, bcf);
      }

      assembler.finalize();

      size_t systsz = assembler.LHS.rows();
      size_t nnz    = assembler.LHS.nonZeros();

      if (verbose())
      {
         std::cout << "Starting linear solver..." << std::endl;
         std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
         std::cout << " * Matrix fill: " << 100.0*double(nnz)/(systsz*systsz) << "%" << std::endl;
      }

      vector_type sol = vector_type::Zero(systsz);

      disk::solvers::pardiso_params<scalar_type> pparams;
      pparams.report_factorization_Mflops = false;
      mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);


      m_postprocess_data.reserve(m_msh.cells_size());

      for (auto& cl : m_msh)
      {
          auto        cb   = make_scalar_monomial_basis(m_msh, cl, m_hdi.cell_degree());
          auto        gr   = make_scalar_hho_laplacian(m_msh, cl, m_hdi);
          auto        stab = make_scalar_hho_stabilization(m_msh, cl, gr.first, m_hdi);
          auto        rhs  = make_rhs(m_msh, cl, cb, lf);
          matrix_type A    = gr.second + stab;

          vector_type locsol = assembler.take_local_data(m_msh, cl, sol, bcf);

          vector_type fullsol = make_scalar_static_decondensation(m_msh, cl, m_hdi, A, rhs, locsol);
          m_postprocess_data.push_back(fullsol);
      }
   }


   // void
   // plot_solution_at_gausspoint(const std::string& filename)
   // {
   //    std::cout << "Compute solution at Gauss points" << std::endl;
   //    visu::Gmesh msh(m_msh.dimension); //creta a mesh

   //    std::vector<visu::Data> data; //create data (not used)
   //    std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point
   //    size_t nb_node =  msh.getNumberofNodes();

   //    size_t cell_i = 0;
   //    for (auto& cl : m_msh)
   //    {
   //       auto x = m_postprocess_data.at(cell_i++);
   //       auto qps = m_bqd.cell_quadrature.integrate(m_msh, cl);
   //       for (auto& qp : qps)
   //       {
   //          auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, qp.point());

   //          scalar_type pot = 0.0;
   //          for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_cell_degree).size(); i++)
   //          pot += phi[i] * x(i);

   //          nb_node += 1;
   //          visu::Node snode = visu::convertPoint(qp.point(), nb_node); //create a node at gauss point
   //          std::vector<double> value = {double(pot)}; // save the solution at gauss point
   //          visu::SubData sdata(value, snode);
   //          subdata.push_back(sdata); // add subdata
   //       }
   //    }

   //    visu::NodeData nodedata(1, 0.0, "sol_scalar", data, subdata); // create and init a nodedata view

   //    nodedata.saveNodeData(filename, msh); // save the view
   // }

   // void
   // plot_conforme_solution(const std::string& filename)
   // {
   //    visu::Gmesh gmsh = visu::convertMesh(m_msh);
   //    auto storage = m_msh.backend_storage();
   //    size_t nb_nodes(gmsh.getNumberofNodes());

   //    //first(number of data at this node), second(cumulated value)
   //    std::vector<std::pair<size_t, scalar_type > > value(nb_nodes, std::make_pair(0, double(0.0)));

   //    size_t cell_i(0);
   //    for (auto& cl : m_msh)
   //    {
   //       vector_dynamic x = m_postprocess_data.at(cell_i++);
   //       auto cell_nodes = visu::cell_nodes(m_msh, cl);
   //       std::vector<visu::Node> new_nodes;
   //       for (size_t i = 0; i < cell_nodes.size(); i++)
   //       {
   //          nb_nodes++;
   //          auto point_ids = cell_nodes[i];
   //          auto pt = storage->points[point_ids];

   //          auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt);

   //          scalar_type pot(0.0);

   //          //compute solution at the node
   //          for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_cell_degree).size(); i++)
   //          pot += phi[i] * x(i);

   //          value[point_ids].first +=1;
   //          value[point_ids].second += pot;
   //       }
   //    }

   //    std::vector<visu::Data> data; //create data
   //    std::vector<visu::SubData> subdata; //create subdata
   //    data.reserve(nb_nodes); // data has a size of nb_node

   //    for(size_t  i_node = 0; i_node < value.size(); i_node++){
   //       std::vector<double> tmp_value(1,value[i_node].second/ double(value[i_node].first));
   //       visu::Data tmp_data(i_node + 1, tmp_value);
   //       data.push_back(tmp_data); //add data
   //    }

   //    visu::NodeData nodedata(1, 0.0, "sol_node", data, subdata); // create and init a nodedata view

   //    nodedata.saveNodeData(filename, gmsh); // save the view
   // }

   // void
   // plot_discontinuous_solution(const std::string& filename)
   // {
   //    const size_t dim = m_msh.dimension;
   //    visu::Gmesh gmsh(dim);
   //    auto storage = m_msh.backend_storage();

   //    std::vector<visu::Data> data; //create data (not used)
   //    std::vector<visu::SubData> subdata; //create subdata to save soution at gauss point

   //    size_t cell_i(0);
   //    size_t nb_nodes(0);
   //    for (auto& cl : m_msh)
   //    {
   //       vector_dynamic x = m_postprocess_data.at(cell_i++);
   //       auto cell_nodes = visu::cell_nodes(m_msh, cl);
   //       std::vector<visu::Node> new_nodes;
   //       for (size_t i = 0; i < cell_nodes.size(); i++)
   //       {
   //          nb_nodes++;
   //          auto point_ids = cell_nodes[i];
   //          auto pt = storage->points[point_ids];

   //          auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt);

   //          std::array<double, 3> coor = {double{0.0}, double{0.0}, double{0.0}};

   //          visu::init_coor(pt, coor);
   //          visu::Node tmp_node(coor, nb_nodes, 0);
   //          new_nodes.push_back(tmp_node);
   //          gmsh.addNode(tmp_node);

   //          // plot magnitude at node
   //          scalar_type pot(0.0);
   //          for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_cell_degree).size(); i++) //compute solution at the node
   //          pot += phi[i] * x(i);

   //          std::vector<double> value = {double(pot)};
   //          visu::Data datatmp(nb_nodes, value);
   //          data.push_back(datatmp);
   //       }
   //       // add new element
   //       visu::add_element(gmsh, new_nodes);
   //    }

   //    visu::NodeData nodedata(1, 0.0, "sol_node", data, subdata); // create and init a nodedata view

   //    nodedata.saveNodeData(filename, gmsh); // save the view
   // }

   // void
   // plot_deformed_conforme(const std::string& filename)
   // {
   //    const size_t DIM = m_msh.dimension;
   //    if(DIM >= 3)
   //    std::cout << "Compute deformed only in 1D or 2D" << '\n';
   //    else {
   //       visu::Gmesh gmsh = visu::convertMesh(m_msh);
   //       auto storage = m_msh.backend_storage();
   //       size_t nb_nodes(gmsh.getNumberofNodes());

   //       //first(number of data at this node), second(cumulated value)
   //       std::vector<std::pair<size_t, scalar_type > > value(nb_nodes, std::make_pair(0, double(0.0)));

   //       size_t cell_i(0);
   //       for (auto& cl : m_msh)
   //       {
   //          vector_dynamic x = m_postprocess_data.at(cell_i++);
   //          auto cell_nodes = visu::cell_nodes(m_msh, cl);
   //          std::vector<visu::Node> new_nodes;
   //          for (size_t i = 0; i < cell_nodes.size(); i++)
   //          {
   //             nb_nodes++;
   //             auto point_ids = cell_nodes[i];
   //             auto pt = storage->points[point_ids];

   //             auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, pt);

   //             scalar_type pot(0.0);

   //             //compute solution at the node
   //             for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_cell_degree).size(); i++)
   //             pot += phi[i] * x(i);

   //             value[point_ids].first +=1;
   //             value[point_ids].second += pot;
   //          }
   //       }

   //       std::vector<visu::Data> data; //create data
   //       std::vector<visu::SubData> subdata; //create subdata
   //       data.reserve(nb_nodes); // data has a size of nb_node

   //       for(size_t  i_node = 0; i_node < value.size(); i_node++){
   //          std::vector<double> tmp_value(3, 0.0);
   //          tmp_value[DIM] = value[i_node].second/ double(value[i_node].first);
   //          visu::Data tmp_data(i_node + 1, tmp_value);
   //          data.push_back(tmp_data); //add data
   //       }

   //       visu::NodeData nodedata(3, 0.0, "sol_node", data, subdata); // create and init a nodedata view

   //       nodedata.saveNodeData(filename, gmsh); // save the view
   //    }
   // }

   void
   plot_deformed_discontinuous(const std::string& filename)
   {
      const size_t DIM = m_msh.dimension;
      if(DIM >= 3)
      std::cout << "Compute deformed only in 1D or 2D" << '\n';
      else {
         gmsh::Gmesh gmsh(DIM);
         auto storage = m_msh.backend_storage();

         std::vector<gmsh::Data>    data;    // create data (not used)
         std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point

         size_t cell_i = 0;
         size_t nb_nodes = 0;
         for (auto& cl : m_msh)
         {
             vector_type          x          = m_postprocess_data.at(cell_i);
             const auto              cell_nodes = disk::cell_nodes(m_msh, cl);
             auto                    cbas       = disk::make_scalar_monomial_basis(m_msh, cl, m_hdi.cell_degree());

             std::vector<gmsh::Node> new_nodes;
             for (size_t i = 0; i < cell_nodes.size(); i++)
             {
                 nb_nodes++;
                 auto point_ids = cell_nodes[i];
                 auto pt        = storage->points[point_ids];

                 const auto phi = cbas.eval_functions(pt);

                 const std::array<double, 3> coor = disk::init_coor(pt);

                 gmsh::Node tmp_node(coor, nb_nodes, 0);
                 new_nodes.push_back(tmp_node);
                 gmsh.addNode(tmp_node);

                 // plot magnitude at node
                 const auto pot = disk::eval(x, phi);

                 std::vector<double> value(3, 0.0);
                 value[DIM + 1 - 1] = {double(pot)};
                 gmsh::Data datatmp(nb_nodes, value);

                 data.push_back(datatmp);
            }
            // add new element
            disk::add_element(gmsh, new_nodes);

            cell_i++;
         }

         gmsh::NodeData nodedata(3, 0.0, "sol_node_deformed", data, subdata); // create and init a nodedata view

         nodedata.saveNodeData(filename, gmsh); // save the view
      }
   }

   void
   saveMesh(const std::string& filename)
   {
       gmsh::Gmesh gmsh = disk::convertMesh(m_msh);
       gmsh.writeGmesh(filename, 2);
   }
};
