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
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "colormanip.h"

#include "../../config.h"

#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "loaders/loader.hpp"
#include "hho/hho.hpp"
#include "ElasticityParameters.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "hyperelasticity_solver.hpp"

struct error_type
{
   size_t  degree;
   size_t nb_dof;
   double h;
   double error_depl;
   double error_grad;
};


struct run_params
{
   size_t  degree;
   int     l;
   bool    verbose;
   size_t n_time_step;
   size_t mindeg;
   size_t maxdeg;
};


template<template<typename, size_t , typename> class Mesh,
typename T, typename Storage>
error_type
run_hyperelasticity_solver(const Mesh<T, 2, Storage>& msh, const run_params& rp, const ElasticityParameters& elas_param)
{
   typedef Mesh<T, 2, Storage> mesh_type;
   typedef static_vector<T, 2> result_type;
   typedef static_matrix<T, 2, 2> result_grad_type;

   T alpha = 0.2;

   auto load = [elas_param, alpha](const point<T,2>& p) -> result_type {
      T lambda = elas_param.lambda;
      T mu = elas_param.mu;

      T num = -lambda * log(9.0 * std::pow(p.x(), 2.0) * std::pow(p.y(), 2.0) + 3.0 * (std::pow(p.x(), 2.0) + std::pow(p.y(), 2.0)) +1.0)
      + lambda + mu;

      T dem1 = 9.0 * std::pow(p.x(), 4.0) + 6.0 * std::pow(p.x(), 2.0) + 1.0;
      T dem2 = 9.0 * std::pow(p.y(), 4.0) + 6.0 * std::pow(p.y(), 2.0) + 1.0;

      T fx = 0.0;//- 6.0 * p.x() * ( num + mu*dem1)/dem1 ;
      T fy = 0.0;//- 6.0 * p.y() * ( num + mu*dem2)/dem2 ;
      return result_type{fx,fy};
   };

   auto solution = [elas_param, alpha](const point<T,2>& p) -> result_type {
      T lambda = elas_param.lambda;
      T fx = (1.0/lambda + alpha) * p.x();
      T fy = (1.0 - alpha/(1.0 + alpha)) * p.y() + /* f(x)= */ alpha * std::pow(p.x(), 1.0);

      return result_type{fx,fy};
   };

   auto gradient = [elas_param, alpha](const point<T,2>& p) -> result_grad_type {
      T lambda = elas_param.lambda;
      result_grad_type grad = result_grad_type::Zero();

      grad(0,0) = (1.0/lambda + alpha);
      grad(1,1) = (1.0 - alpha/(1.0 + alpha));
      grad(1,0) = /* f'(x)= */ alpha;

      return grad;
   };

   auto neumann = [elas_param](const point<T,2>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.0;

      return result_type{fx,fy};
   };

   std::vector<size_t> boundary_neumann(0);

   hyperelasticity_solver<Mesh, T, 2, Storage,  point<T, 2> > nl(msh, rp.degree, elas_param);
   nl.verbose(false);


   nl.compute_initial_state(boundary_neumann);

   solve_info solve_info = nl.compute(load, solution, neumann, boundary_neumann, rp.n_time_step);

   error_type error;
   error.h = average_diameter(msh);
   error.degree = rp.degree;
   error.nb_dof = nl.getDofs();
   error.error_depl = 10E6;
   error.error_grad = 10E6;

   if(nl.test_convergence()){
      error.error_depl = nl.compute_l2_error(solution);
      error.error_grad = nl.compute_l2_gradient_error(gradient);
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
      std::cout << "-----------------------------------------------------------------------------------" << std::endl;
      std::cout << "| Size mesh  | Displacement | Convergence |  Gradient  | Convergence |    Total   |" << std::endl;
      std::cout << "|    h       |   L2 error   |     rate    |  L2 error  |     rate    | faces DOF  |" << std::endl;
      std::cout << "-----------------------------------------------------------------------------------" << std::endl;


      std::string s_dof = " " + std::to_string(error[0].nb_dof) + "                  ";
      s_dof.resize(10);

      std::cout << "| " <<  error[0].h << " |  " << error[0].error_depl << "  | " << "     -     " << " | " <<
      error[0].error_grad <<  " | " << "     -     "  <<  " | " << s_dof  <<  " |" << std::endl;

      for(size_t i = 1; i < error.size(); i++){
         s_dof = " " + std::to_string(error[i].nb_dof) + "                  ";
         s_dof.resize(10);
         double rate_depl = (log10(error[i-1].error_depl) - log10(error[i].error_depl))/(log10(error[i-1].h) - log10(error[i].h));
         double rate_grad = (log10(error[i-1].error_grad) - log10(error[i].error_grad))/(log10(error[i-1].h) - log10(error[i].h));

         std::cout << "| " <<  error[i].h << " |  " << error[i].error_depl << "  |  " << rate_depl << " | " <<
         error[i].error_grad <<  " |  " << rate_grad  <<  " | " << s_dof <<  " |" << std::endl;
      }

      std::cout << "-----------------------------------------------------------------------------------" << std::endl;
      std::cout << "  " <<std::endl;
      std::cout.flags( f );
   }
   else
      std::cout << "The file error is empty" << std::endl;
}

template< typename T>
void test_triangles_fvca5(const run_params& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 5;

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
void test_triangles_netgen(const run_params& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 5;

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
void test_hexagons(const run_params& rp, const ElasticityParameters& elas_param)
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
void test_kershaws(const run_params& rp, const ElasticityParameters& elas_param)
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
void test_quads_fvca5(const run_params& rp, const ElasticityParameters& elas_param)
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
void test_quads_diskpp(const run_params& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/2D_quads/diskpp/testmesh-2-2.quad");
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



int main(int argc, char **argv)
{
   using RealType = double;

   char    *mesh_filename  = nullptr;
   char    *plot_filename  = nullptr;
   int     degree          = 1;
   int     l               = 0;
   int     elems_1d        = 8;

   run_params rp;
   rp.degree   = 1;
   rp.l        = 0;
   rp.verbose  = false;
   rp.n_time_step = 1;

   ElasticityParameters param = ElasticityParameters();
   param.lambda = 100.0;
   param.mu = 1.0;
   param.tau = 1000.0;
   param.adaptative_stab = false;
   param.type_law = 1;

   int ch;

   while ( (ch = getopt(argc, argv, "k:l:n:p:v:t")) != -1 )
   {
      switch(ch)
      {
         case 'k':
            degree = atoi(optarg);
            if (degree < 0)
            {
               std::cout << "Degree must be positive. Falling back to 1." << std::endl;
               degree = 1;
            }
            rp.degree = degree;
            break;

         case 'l':
            rp.l = atoi(optarg);
            if (l < -1 or l > 1)
            {
               std::cout << "l can be -1, 0 or 1. Falling back to 0." << std::endl;
               rp.l = 0;
            }
            break;

         case 'n':
            rp.n_time_step = atoi(optarg);
            if (rp.n_time_step == 0)
            {
               std::cout << "Number of time step must be positive. Falling back to 1." << std::endl;
               rp.n_time_step = 1;
            }
            break;


         case 't':
            param.tau = atof(optarg);
            break;

         case 'v':
            rp.verbose = true;
            break;

         case 'h':
         case '?':
         default:
            std::cout << "wrong arguments" << std::endl;
            exit(1);
      }
   }

   argc -= optind;
   argv += optind;


   timecounter tc;

   std::cout << " Test convergence rates in 2D for: "<< std::endl;
   std::cout << " ** Degree k = " << rp.degree << std::endl;
   std::cout << " ** Stab tau = " << param.tau << std::endl;
   std::cout << " ** mu = " << param.mu << std::endl;
   std::cout << " ** lambda = " << param.lambda << std::endl;
   std::cout << " "<< std::endl;

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
