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


#include "../exemple_visualisation/visualisation/gmshDisk.hpp"
#include "../exemple_visualisation/visualisation/gmshConvertMesh.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>


struct assembly_info
{
   size_t  linear_system_size;
   double  time_gradrec, time_statcond, time_stab, time_assembly;
};

struct solver_info
{
   double  time_solver;
};

struct postprocess_info
{
   double  time_postprocess;
};

struct LaplacianParameters {
   LaplacianParameters(double _lambda = 1.)
   : lambda(_lambda)
   {
      // Do nothing
   }
   
   double lambda;
};


template<typename Mesh>
class vector_laplacian_solver
{
   typedef Mesh                     mesh_type;
   typedef typename mesh_type::scalar_type            scalar_type;
   typedef typename mesh_type::cell                   cell_type;
   typedef typename mesh_type::face                   face_type;
   
   typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
   typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;
   typedef disk::quadrature<mesh_type, cell_type>      matrix_quadrature_type;
   
   typedef disk::scaled_monomial_vector_basis<mesh_type, cell_type>    cell_basis_type;
   typedef disk::scaled_monomial_vector_basis<mesh_type, face_type>    face_basis_type;
   typedef disk::scaled_monomial_matrix_basis<mesh_type, cell_type>    matrix_basis_type;
   
   typedef dynamic_matrix<scalar_type>         matrix_dynamic;
   typedef dynamic_vector<scalar_type>         vector_dynamic;
   
   typedef disk::gradient_reconstruction_elas_full< mesh_type,
                                                cell_basis_type, cell_quadrature_type,
                                                face_basis_type, face_quadrature_type,
                                                matrix_basis_type, matrix_quadrature_type>          gradrec_type;
   
   typedef disk::elas_like_stabilization_l2<   mesh_type,
                                             cell_basis_type, cell_quadrature_type,
                                             face_basis_type, face_quadrature_type>                stab_type;
   
   typedef disk::diffusion_like_static_condensation<  mesh_type,
                                                   cell_basis_type, cell_quadrature_type,
                                                   face_basis_type, face_quadrature_type>    statcond_type;
   
   typedef disk::assembler_elas<mesh_type, face_basis_type, face_quadrature_type>               assembler_type;
   
   typedef disk::projector_elas<mesh_type, cell_basis_type, cell_quadrature_type,
                               face_basis_type, face_quadrature_type>                projector_type;
   
   typename assembler_type::sparse_matrix_type     m_system_matrix;
   typename assembler_type::vector_type            m_system_rhs, m_system_solution;
   
   
   size_t m_cell_degree, m_face_degree, m_degree;
   
   const mesh_type& m_msh;
   
   std::vector<vector_dynamic>                    m_solution_data;
   
   bool m_verbose;
   
   LaplacianParameters m_laplacian_parameters;
   
public:
   vector_laplacian_solver(const mesh_type& msh, size_t degree, int l = 0)
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
      
      m_cell_degree = degree + l;
      m_face_degree = degree;
      m_degree = degree;
      
      m_laplacian_parameters.lambda = 1.0;
      
   }
   
   vector_laplacian_solver(const mesh_type& msh, size_t degree, const LaplacianParameters data, int l = 0)
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
      
      m_cell_degree = degree + l;
      m_face_degree = degree;
      m_degree = degree;
      
      m_laplacian_parameters.lambda = data.lambda;
      
   }
   
   void
   changeLaplacianParameters(const LaplacianParameters data)
   {
      m_laplacian_parameters.lambda = data.lambda;
   }
   
   bool    verbose(void) const     { return m_verbose; }
   void    verbose(bool v)         { m_verbose = v; }
   
   
   template<typename LoadFunction, typename BoundaryConditionFunction>
   assembly_info
   assemble(const LoadFunction& lf, const BoundaryConditionFunction& bcf)
   {
      auto gradrec     = gradrec_type(m_degree);
      auto stab         = stab_type(m_degree);
      auto statcond     = statcond_type(m_degree);
      auto assembler    = assembler_type(m_msh, m_face_degree);
      
      assembly_info ai;
      bzero(&ai, sizeof(ai));
      
      timecounter tc;
      
      for (auto& cl : m_msh)
      {
         tc.tic();
         gradrec.compute(m_msh, cl);
         tc.toc();
         ai.time_gradrec += tc.to_double();
         
         tc.tic();
         stab.compute(m_msh, cl, gradrec.oper);
         tc.toc();
         ai.time_stab += tc.to_double();
         
         assert(gradrec.data.rows() == stab.data.rows());
         assert(gradrec.data.cols() == stab.data.cols());
         
         tc.tic();
         auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(m_msh, cl, lf, m_cell_degree);
         dynamic_matrix<scalar_type> loc = m_laplacian_parameters.lambda * (gradrec.data + stab.data);
         auto scnp = statcond.compute(m_msh, cl, loc, cell_rhs);
         tc.toc();
         ai.time_statcond += tc.to_double();
         
         assembler.assemble(m_msh, cl, scnp);
      }
      
      assembler.impose_boundary_conditions(m_msh, bcf);
      assembler.finalize(m_system_matrix, m_system_rhs);
      
      ai.linear_system_size = m_system_matrix.rows();
      ai.time_assembly = ai.time_gradrec + ai.time_stab + ai.time_statcond;
      return ai;
   }
   
   solver_info
   solve(void)
   {
      #ifdef HAVE_INTEL_MKL
      Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
      #else
      Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
      #endif
      
      solver_info si;
      
      size_t systsz = m_system_matrix.rows();
      size_t nnz = m_system_matrix.nonZeros();
      
      if (verbose())
      {
         std::cout << "Starting linear solver..." << std::endl;
         std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
         std::cout << " * Matrix fill: " << 100.0*double(nnz)/(systsz*systsz) << "%" << std::endl;
      }
      
      //       if(m_verbose){
      //          std::cout << "** Numbers of dofs: " << total_dof + total_lagr  << std::endl;
      //          std::cout << "** including " << total_lagr << " Lagrange multipliers"  << std::endl;
      //          std::cout << "** After static condensation: "  << std::endl;
      //          std::cout << "** Numbers of dofs: " << m_msh.faces_size() * num_face_dofs + total_lagr  << std::endl;
      //          std::cout << "** including " << total_lagr << " Lagrange multipliers"  << std::endl;
      //       }
      
      timecounter tc;
      
      tc.tic();
      solver.analyzePattern(m_system_matrix);
      solver.factorize(m_system_matrix);
      m_system_solution = solver.solve(m_system_rhs);
      tc.toc();
      si.time_solver = tc.to_double();
      
      return si;
   }
   
   template<typename LoadFunction>
   postprocess_info
   postprocess(const LoadFunction& lf)
   {
      auto gradrec      = gradrec_type(m_degree);
      auto stab         = stab_type(m_degree);
      auto statcond     = statcond_type(m_degree);
      
      face_basis_type face_basis(m_face_degree);
      size_t fbs = face_basis.size();
      
      postprocess_info pi;
      
      m_solution_data.reserve(m_msh.cells_size());
      
      timecounter tc;
      tc.tic();
      for (auto& cl : m_msh)
      {
         auto fcs = faces(m_msh, cl);
         auto num_faces = fcs.size();
         
         dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(num_faces*fbs);
         
         for (size_t face_i = 0; face_i < num_faces; face_i++)
         {
            auto fc = fcs[face_i];
            auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), fc);
            if (!eid.first)
               throw std::invalid_argument("This is a bug: face not found");
            
            auto face_id = eid.second;
            
            dynamic_vector<scalar_type> xF = dynamic_vector<scalar_type>::Zero(fbs);
            xF = m_system_solution.block(face_id * fbs, 0, fbs, 1);
            xFs.block(face_i * fbs, 0, fbs, 1) = xF;
         }
         
         gradrec.compute(m_msh, cl);
         stab.compute(m_msh, cl, gradrec.oper);
         dynamic_matrix<scalar_type> loc = m_laplacian_parameters.lambda * (gradrec.data + stab.data);
         auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(m_msh, cl, lf, m_cell_degree);
         dynamic_vector<scalar_type> x = statcond.recover(m_msh, cl, loc, cell_rhs, xFs);
         m_solution_data.push_back(x);
      }
      tc.toc();
      
      pi.time_postprocess = tc.to_double();
      
      return pi;
   }
   
   
   
   
   template<typename AnalyticalSolution>
   scalar_type
   compute_l2_error(const AnalyticalSolution& as)
   {
      scalar_type err_dof = scalar_type{0.0};
      
      projector_type projk(m_degree);
      
      size_t i = 0;
      
      for (auto& cl : m_msh)
      {
         auto x = m_solution_data.at(i++);
         dynamic_vector<scalar_type> true_dof = projk.compute_cell(m_msh, cl, as);
         dynamic_vector<scalar_type> comp_dof = x.block(0,0,true_dof.size(), 1);
         dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
         err_dof += diff_dof.dot(projk.cell_mm * diff_dof);
      }
      
      return sqrt(err_dof);
   }
   
   /*
    v o*id
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
saveMesh(const std::string& filename)
{
visu::Gmesh gmsh = visu::convertMesh(m_msh);
gmsh.writeGmesh(filename, DIM);
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
}*/
   
};