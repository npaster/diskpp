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
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "colormanip.h"

#include "config.h"

#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "loaders/loader.hpp"
#include "hho/hho.hpp"
#include "ElasticityParameters.hpp"
#include "mechanics/BoundaryConditions.hpp"
#include "Parameters.hpp"
#include "Informations.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "hyperelasticity_solver.hpp"

/* TEST CONVERGENCE RATES FOR HYPERELASTICITY PROBLEMS
* c.f. Article: Hybrid High-Order methods for finite deformations of hyperlastic problems
* M. Abbas, A. Ern and N. Pignet
* Computational Mechanics
* DOI:
*/

struct error_type
{
   size_t  degree;
   size_t nb_dof;
   double h;
   double error_depl;
   double error_grad;
   double error_energy;
   double error_PK1;
};


void
usage(const char *progname)
{
   printf("Usage: %s <options> <filename>\n\n", progname);
   printf("    -2: test 2D mesh (default)\n");
   printf("    -3: test 3D mesh\n");
   printf("    -k: face degree (>0)\n");
   printf("    -l: difference beetween cell and face degree (-1 <= l <= 1) \n");
   printf("    -n: number of time step (>0)\n");
   printf("    -m: number of sublevel time step (>0)\n");
   printf("    -v: verbose\n");
   printf("    -t: stabilization parameter (>0)\n");
   printf("    -s: adaptative stabilisation\n");
}


template<template<typename, size_t , typename> class Mesh,
typename T, typename Storage>
error_type
run_hyperelasticity_solver(const Mesh<T, 2, Storage>& msh, const ParamRun<T>& rp, const ElasticityParameters& elas_param)
{
   typedef static_vector<T, 2> result_type;
   typedef static_matrix<T, 2, 2> result_grad_type;

   // Load ans solutions
   const T alpha = 0.3;

   auto load = [elas_param, alpha](const point<T,2>& p) -> result_type {
      T mu = elas_param.mu;

      T fx = 0.0;
      T fy = 8*mu * alpha * M_PI* M_PI* cos(2*M_PI*p.x());
      return result_type{fx,fy};
   };

   auto solution = [elas_param, alpha](const point<T,2>& p) -> result_type {
      T lambda = elas_param.lambda;
      T fx = (1.0/lambda + alpha) * p.x();
      T fy = (1.0/lambda - alpha/(1.0 + alpha)) * p.y() + /* f(x)= */ 2*alpha * (cos(2*M_PI*p.x()) -1.0);

      return result_type{fx,fy};
   };

   auto gradient = [elas_param, alpha](const point<T,2>& p) -> result_grad_type {
      T lambda = elas_param.lambda;
      result_grad_type grad = result_grad_type::Zero();

      grad(0,0) = (1.0/lambda + alpha);
      grad(1,1) = (1.0/lambda - alpha/(1.0 + alpha));
      grad(1,0) = /* f'(x)= */ -4*alpha * M_PI* sin(2*M_PI*p.x());

      return grad;
   };

   auto neumann = [elas_param](const point<T,2>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.0;

      return result_type{fx,fy};
   };

   const std::vector<BoundaryType> boundary_neumann = {};
   const std::vector<BoundaryType> boundary_dirichlet = {};

   // Solve the problem
   hyperelasticity_solver<Mesh, T, 2, Storage,  point<T, 2> >
      nl(msh, rp, elas_param, boundary_neumann, boundary_dirichlet);

   nl.compute_initial_state();

   const SolverInfo solve_info = nl.compute(load, solution, neumann);

   // Compute errors
   error_type error;
   error.h = average_diameter(msh);
   error.degree = rp.m_face_degree;
   error.nb_dof = nl.getDofs();
   error.error_depl = 10E6;
   error.error_grad = 10E6;
   error.error_energy = 10E6;
   error.error_PK1 = 10E6;

   if(nl.test_convergence()){
      error.error_depl = nl.compute_l2_error_displacement(solution);
      try{
         error.error_grad = nl.compute_l2_error_gradient(gradient);
      }
      catch(const std::invalid_argument& ia){
         std::cerr << "Invalid argument: " << ia.what() << " in compute_l2_error_gradient" << std::endl;
      }
      try{
         error.error_energy = nl.compute_l2_error_energy(gradient);
      }
      catch(const std::invalid_argument& ia){
         std::cerr << "Invalid argument: " << ia.what() << " in compute_l2_error_energy" << std::endl;
      }
      try{
         error.error_PK1 = nl.compute_l2_error_PK1(gradient);
      }
      catch(const std::invalid_argument& ia){
         std::cerr << "Invalid argument: " << ia.what() << " in compute_l2_error_PK1" << std::endl;
      }
   }

   if(nl.verbose()){
      std::cout << " " << std::endl;
      std::cout << "------------------------------------------------------- " << std::endl;
      std::cout << "Summaring: " << std::endl;
      std::cout << "Total Newton's iterations: " << solve_info.m_iter << " in " << solve_info.m_time_step << " load increments" << std::endl;
      std::cout << "Total time to solve the problem: " << solve_info.m_time_solver << " sec" << std::endl;
      std::cout << "**** Assembly time: " << solve_info.m_newton_info.m_assembly_info.m_time_assembly << " sec" << std::endl;
      std::cout << "****** Gradient reconstruction: " << solve_info.m_newton_info.m_assembly_info.m_time_gradrec << " sec" << std::endl;
      std::cout << "****** Stabilisation: " << solve_info.m_newton_info.m_assembly_info.m_time_stab << " sec" << std::endl;
      std::cout << "****** Elementary computation: " << solve_info.m_newton_info.m_assembly_info.m_time_elem << " sec" << std::endl;
      std::cout << "       *** Behavior computation: " << solve_info.m_newton_info.m_assembly_info.m_time_law << " sec" << std::endl;
      std::cout << "****** Static condensation: " << solve_info.m_newton_info.m_assembly_info.m_time_statcond << " sec" << std::endl;
      std::cout << "**** Postprocess time: " << solve_info.m_newton_info.m_assembly_info.m_time_postpro << " sec" << std::endl;
      std::cout << "**** Solver time: " << solve_info.m_newton_info.m_solve_info.m_time_solve << " sec" << std::endl;
      std::cout << "------------------------------------------------------- " << std::endl;
      std::cout << " " << std::endl;
   }

   return error;
}


template<template<typename, size_t, typename> class Mesh,
typename T, typename Storage>
error_type
run_hyperelasticity_solver(const Mesh<T, 3, Storage>& msh, const ParamRun<T>& rp, const ElasticityParameters& elas_param)
{
   typedef static_vector<T, 3> result_type;
   typedef static_matrix<T, 3, 3> result_grad_type;

   // Load and solutions
   const T alpha = 0.2;
   const T beta = 0.2;
   const T factor = 0.5;

   auto load = [elas_param, alpha,beta, factor](const point<T,3>& p) -> result_type {
      T fx = M_PI*M_PI*elas_param.mu*alpha * sin(M_PI*p.y());
      T fy = 0.0;
      T fz = M_PI*M_PI*elas_param.mu*beta * sin(M_PI*p.x());
      return factor*result_type{fx, fy, fz};
   };

   auto solution = [elas_param, alpha, beta, factor](const point<T,3>& p) -> result_type {
      T lambda = elas_param.lambda;
      T gamma = alpha + beta + alpha * beta;
      T fx = (1.0/lambda + alpha) * p.x() + /* f(Y)= */ alpha * sin(M_PI*p.y());
      T fy = (1.0/lambda - gamma/(1.0 + gamma)) * p.y();
      T fz = (1.0/lambda + beta) * p.z() + /* g(X)= */ beta * sin(M_PI*p.x()) + /* h(Y)= */ 0.0;

      return factor*result_type{fx,fy,fz};
   };

   auto gradient = [elas_param, alpha, beta, factor](const point<T,3>& p) -> result_grad_type {
      T lambda = elas_param.lambda;
      T gamma = alpha + beta + alpha * beta;
      result_grad_type grad = result_grad_type::Zero();

      grad(0,0) = (1.0/lambda + alpha);
      grad(0,1) = /* f'(Y)= */ M_PI*alpha * cos(M_PI*p.y());

      grad(1,1) = (1.0/lambda - gamma/(1.0 + gamma));

      grad(2,0) = /* g'(X)= */ M_PI*beta* cos(M_PI*p.x());
      grad(2,1) = /* h'(Y)= */ 0.0;
      grad(2,2) = (1.0/lambda + beta);

      return factor*grad;
   };


   auto neumann = [elas_param](const point<T,3>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.0;
      T fz = 0.0;

      return result_type{fx,fy,fz};
   };

   const std::vector<BoundaryType> boundary_neumann = {};
   const std::vector<BoundaryType> boundary_dirichlet = {};

   // Solve
   hyperelasticity_solver<Mesh, T, 3, Storage,  point<T, 3> >
      nl(msh, rp, elas_param, boundary_neumann, boundary_dirichlet);

   nl.compute_initial_state();

   const SolverInfo solve_info = nl.compute(load, solution, neumann);

   // Compute errors
   error_type error;
   error.h = average_diameter(msh);
   error.degree = rp.m_face_degree;
   error.nb_dof = nl.getDofs();
   error.error_depl = 10E6;
   error.error_grad = 10E6;
   error.error_energy = 10E6;
   error.error_PK1 = 10E6;

   if(nl.test_convergence()){
      error.error_depl = nl.compute_l2_error_displacement(solution);
      try{
         error.error_grad = nl.compute_l2_error_gradient(gradient);
      }
      catch(const std::invalid_argument& ia){
         std::cerr << "Invalid argument: " << ia.what() << " in compute_l2_error_gradient" << std::endl;
      }
      try{
         error.error_energy = nl.compute_l2_error_energy(gradient);
      }
      catch(const std::invalid_argument& ia){
         std::cerr << "Invalid argument: " << ia.what() << " in compute_l2_error_energy" << std::endl;
      }
      try{
         error.error_PK1 = nl.compute_l2_error_PK1(gradient);
      }
      catch(const std::invalid_argument& ia){
         std::cerr << "Invalid argument: " << ia.what() << " in compute_l2_error_PK1" << std::endl;
      }
      nl.test_pt(solution, gradient);
   }

   if(nl.verbose()){
      std::cout << " " << std::endl;
      std::cout << "------------------------------------------------------- " << std::endl;
      std::cout << "Summaring: " << std::endl;
      std::cout << "Total Newton's iterations: " << solve_info.m_iter << " in " << solve_info.m_time_step << " load increments" << std::endl;
      std::cout << "Total time to solve the problem: " << solve_info.m_time_solver << " sec" << std::endl;
      std::cout << "**** Assembly time: " << solve_info.m_newton_info.m_assembly_info.m_time_assembly << " sec" << std::endl;
      std::cout << "****** Gradient reconstruction: " << solve_info.m_newton_info.m_assembly_info.m_time_gradrec << " sec" << std::endl;
      std::cout << "****** Stabilisation: " << solve_info.m_newton_info.m_assembly_info.m_time_stab << " sec" << std::endl;
      std::cout << "****** Elementary computation: " << solve_info.m_newton_info.m_assembly_info.m_time_elem << " sec" << std::endl;
      std::cout << "       *** Behavior computation: " << solve_info.m_newton_info.m_assembly_info.m_time_law << " sec" << std::endl;
      std::cout << "****** Static condensation: " << solve_info.m_newton_info.m_assembly_info.m_time_statcond << " sec" << std::endl;
      std::cout << "**** Postprocess time: " << solve_info.m_newton_info.m_assembly_info.m_time_postpro << " sec" << std::endl;
      std::cout << "**** Solver time: " << solve_info.m_newton_info.m_solve_info.m_time_solve << " sec" << std::endl;
      std::cout << "------------------------------------------------------- " << std::endl;
      std::cout << " " << std::endl;
   }

   return error;
}

void
printResults(const std::vector<error_type>& error)
{
   if(error.size() > 0){
      std::ios::fmtflags f( std::cout.flags() );
      std::cout.precision(4);
      std::cout.setf(std::iostream::scientific, std::iostream::floatfield);

      std::cout << "Convergence test for k = " << error[0].degree << std::endl;
      std::cout << "-----------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
      std::cout << "| Size mesh  | Displacement | Convergence |  Gradient  | Convergence |  Energy    | Convergence |     PK1    | Convergence |    Total   |" << std::endl;
      std::cout << "|    h       |   L2 error   |     rate    |  L2 error  |     rate    |  L2 error  |     rate    |  L2 error  |     rate    | faces DOF  |" << std::endl;
      std::cout << "-----------------------------------------------------------------------------------------------------------------------------------------" << std::endl;


      std::string s_dof = " " + std::to_string(error[0].nb_dof) + "                  ";
      s_dof.resize(10);

      std::cout << "| " <<  error[0].h << " |  " << error[0].error_depl << "  | " << "     -     " << " | " <<
      error[0].error_grad <<  " | " << "     -     " << " | " <<
      error[0].error_energy <<  " | " << "     -     "  <<  " | " <<
      error[0].error_PK1 <<  " | " << "     -     "  <<  " | "<< s_dof  <<  " |" << std::endl;

      for(size_t i = 1; i < error.size(); i++){
         s_dof = " " + std::to_string(error[i].nb_dof) + "                  ";
         s_dof.resize(10);
         const double rate_depl = (log10(error[i-1].error_depl) - log10(error[i].error_depl))/(log10(error[i-1].h) - log10(error[i].h));
         const double rate_grad = (log10(error[i-1].error_grad) - log10(error[i].error_grad))/(log10(error[i-1].h) - log10(error[i].h));
         const double rate_ener = (log10(error[i-1].error_energy) - log10(error[i].error_energy))/(log10(error[i-1].h) - log10(error[i].h));
         const double rate_PK1 = (log10(error[i-1].error_PK1) - log10(error[i].error_PK1))/(log10(error[i-1].h) - log10(error[i].h));

         std::cout << "| " <<  error[i].h << " |  " << error[i].error_depl << "  |  " << rate_depl << " | " <<
         error[i].error_grad <<  " |  " << rate_grad << " | " <<
         error[i].error_energy <<  " |  " << rate_ener   <<  " | " <<
         error[i].error_PK1 <<  " |  " << rate_PK1   <<  " | "<< s_dof <<  " |" << std::endl;
      }

      std::cout << "-----------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
      std::cout << "  " <<std::endl;
      std::cout.flags( f );
   }
   else
      std::cout << "The file error is empty" << std::endl;
}

// Mesh loaders
template< typename T>
void test_triangles_fvca5(const ParamRun<T>& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_1.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_2.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_3.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_4.typ1");
   paths.push_back("../meshes/2D_triangles/fvca5/mesh1_5.typ1");

   std::vector<error_type> error_sumup;

   for(size_t i = 0; i < runs; i++){
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_triangles_netgen(const ParamRun<T>& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri01.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri02.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri03.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri04.mesh2d");
   paths.push_back("../diskpp/meshes/2D_triangles/netgen/tri05.mesh2d");

   std::vector<error_type> error_sumup;

   for(size_t i = 0; i < runs; i++){
      auto msh = disk::load_netgen_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param));
   }
   printResults(error_sumup);
}



template< typename T>
void test_hexagons(const ParamRun<T>& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_1.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_2.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_3.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_4.typ1");
   paths.push_back("../diskpp/meshes/2D_hex/fvca5/hexagonal_5.typ1");

   std::vector<error_type> error_sumup;

   for(size_t i = 0; i < runs; i++){
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_kershaws(const ParamRun<T>& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_1.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_2.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_3.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_4.typ1");
   paths.push_back("../diskpp/meshes/2D_kershaw/fvca5/mesh4_1_5.typ1");

   std::vector<error_type> error_sumup;

   for(size_t i = 0; i < runs; i++){
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_quads_fvca5(const ParamRun<T>& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_1.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_2.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_3.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_4.typ1");
   paths.push_back("../diskpp/meshes/2D_quads/fvca5/mesh2_5.typ1");

   std::vector<error_type> error_sumup;

   for(size_t i = 0; i < runs; i++){
      auto msh = disk::load_fvca5_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_quads_diskpp(const ParamRun<T>& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-4-4.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-8-8.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-16-16.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-32-32.quad");
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-256-256.quad");

   std::vector<error_type> error_sumup;

   for(size_t i = 0; i < runs; i++){
      auto msh = disk::load_cartesian_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param));
   }
   printResults(error_sumup);
}

template< typename T>
void test_hexahedra_diskpp(const ParamRun<T>& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-2-2-2.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-4-4-4.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-8-8-8.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-16-16-16.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-32-32-32.hex");

   std::vector<error_type> error_sumup;

   for(int i = 0; i < runs; i++){
      auto msh = disk::load_cartesian_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_hexahedra_fvca6(const ParamRun<T>& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_2x2x2.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_4x4x4.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_8x8x8.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_16x16x16.msh");
   paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_32x32x32.hex");

   std::vector<error_type> error_sumup;

   for(int i = 0; i < runs; i++){
      auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param));
   }
   printResults(error_sumup);
}

template< typename T>
void test_tetrahedra_netgen(const ParamRun<T>& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 6;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet0.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/cube336.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet1.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet2.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet3.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet4.mesh");

   std::vector<error_type> error_sumup;

   for(int i = 0; i < runs; i++){
      auto msh = disk::load_netgen_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_polyhedra_fvca6(const ParamRun<T>& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 3;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_10.msh");
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_20.msh");
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_30.msh");
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_40.msh");

   std::vector<error_type> error_sumup;

   for(int i = 0; i < runs; i++){
      auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_tetrahedra_fvca6(const ParamRun<T>& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.1.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.2.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.3.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.4.msh");

   std::vector<error_type> error_sumup;

   for(int i = 0; i < runs; i++){
      auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param));
   }
   printResults(error_sumup);
}

// Main
int main(int argc, char **argv)
{
   using RealType = double;

   int     face_degree          = 1;
   int     cell_degree          = 0;
   int     grad_degree          = 0;
   int     l               = 0;
   int     sublevel        = 1;
   bool three_dimensions = false;

   ParamRun<RealType> rp;
   rp.m_sublevel = 4;
   rp.m_compute_energy = false;

   ElasticityParameters param = ElasticityParameters();
   param.lambda = 10.0;
   param.mu = 1.0;
   param.type_law = 1;

   int ch;

   while ( (ch = getopt(argc, argv, "23b:g:k:l:m:n:r:s:t:v:w")) != -1 )
   {
      switch(ch)
      {
         case '2': three_dimensions = false; break;

         case '3': three_dimensions = true; break;

         case 'g':
            grad_degree = atoi(optarg);
            if (grad_degree <= 0)
            {
                std::cout << "'grad degree' must be positive. Falling back to 1." << std::endl;
                grad_degree = 1;
            }
            rp.m_grad_degree = grad_degree;
            break;

         case 'k':
             face_degree = atoi(optarg);
             if (face_degree <= 0)
             {
                 std::cout << "Degree must be positive. Falling back to 1." << std::endl;
                 face_degree = 1;
             }
             rp.m_face_degree = face_degree;
             if(cell_degree == 0) rp.m_cell_degree = face_degree;
             if(grad_degree == 0) rp.m_grad_degree = face_degree;

             break;

         case 'l':
             l = atoi(optarg);
             if (l < -1 or l > 1)
             {
                 std::cout << "l can be -1, 0 or 1. Falling back to 0." << std::endl;
                 l = 0;
             }
             rp.m_l = l;
             rp.m_cell_degree = rp.m_face_degree + l;
             break;

         case 'm':
            sublevel = atoi(optarg);
            if (sublevel <= 0)
            {
               std::cout << "Number of sublevel time step must be positive. Falling back to 1." << std::endl;
               sublevel = 1;
            }
            rp.m_sublevel = sublevel;
            break;

         case 'r':
            rp.readParameters(optarg); break;


         case 'v': rp.m_verbose = true; break;

         case 'w': rp.m_stab = false; break;

         case '?':
         default:
            std::cout << "wrong arguments" << std::endl;
            usage(argv[0]);
            exit(1);
      }
   }

   argc -= optind;
   argv += optind;

   timecounter tc;

   std::cout << " Test convergence rates for: "<< std::endl;
   std::cout << " ** Face_Degree = " << rp.m_face_degree << std::endl;
   std::cout << " ** Cell_Degree  = " << rp.m_cell_degree << std::endl;
   std::cout << " ** Grad_Degree  = " << rp.m_grad_degree << std::endl;
   std::cout << " ** mu = " << param.mu << std::endl;
   std::cout << " ** lambda = " << param.lambda << std::endl;
   std::cout << " "<< std::endl;

   if(three_dimensions){
      tc.tic();
      std::cout << "-Tetrahedras fvca6:" << std::endl;
      //test_tetrahedra_fvca6<RealType>(rp, param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout <<  "-Tetrahedras netgen:" << std::endl;
      test_tetrahedra_netgen<RealType>(rp, param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout << "-Hexahedras fvca6:"  << std::endl;
      //test_hexahedra_fvca6<RealType>(rp, param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout << "-Hexahedras diskpp:"  << std::endl;
      //test_hexahedra_diskpp<RealType>(rp, param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;


      tc.tic();
      std::cout << "-Polyhedra:"  << std::endl;
      //test_polyhedra_fvca6<RealType>(rp, param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;
   }
   else{

      tc.tic();
      std::cout << "-Triangles fvca5:" << std::endl;
      test_triangles_fvca5<RealType>(rp, param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout <<  "-Triangles netgen:" << std::endl;
      test_triangles_netgen<RealType>(rp, param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout << "-Quadrangles fvca5:"  << std::endl;
      test_quads_fvca5<RealType>(rp, param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

      tc.tic();
      std::cout << "-Quadrangles diskpp:"  << std::endl;
      test_quads_diskpp<RealType>(rp, param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;


      tc.tic();
      std::cout << "-Hexagons:"  << std::endl;
      test_hexagons<RealType>(rp, param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;


      tc.tic();
      std::cout << "-Kershaws:"  << std::endl;
      test_kershaws<RealType>(rp, param);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;

   }

}
