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



template<template<typename, size_t, typename> class Mesh,
typename T, typename Storage>
error_type
run_hyperelasticity_solver(const Mesh<T, 3, Storage>& msh, const run_params& rp, const ElasticityParameters& elas_param)
{
   typedef Mesh<T, 3, Storage> mesh_type;
   typedef static_vector<T, 3> result_type;
   typedef static_matrix<T, 3, 3> result_grad_type;

   // auto load = [](const point<T,3>& p) -> result_type {
   //  !! A calculer
   //    return result_type{0.0,0.0,0.0};
   // };
   //
   // auto solution = [](const point<T,3>& p) -> result_type {
   //    T fx = -1.75 * p.y() * (1-p.y()) * p.z() * (1.0 - p.z()) * cos(M_PI * p.x());
   //    T fy = -1.75 * p.x() * (1-p.x()) * p.z() * (1.0 - p.z()) * cos(M_PI * p.y());
   //    T fz = -0.12 * p.z() * p.z() * (1.0 - cos(2.0 * M_PI * p.x())) * (1.0 - cos(2.0*M_PI * p.y())) + 0.15 * p.z();
   //    return result_type{fx,fy,fz};
   // };
   //
   //
   // auto gradient = [](const point<T,3>& p) -> static_matrix<T, 3, 3> {
   //    static_matrix<T, 3, 3> g;
   //    g(1,1) = -1.75 * M_PI * p.y() * (1-p.y()) * p.z() * (1.0 - p.z()) * sin(M_PI * p.x());
   //    g(1,2) = -1.75 *(1-2 * p.y()) * p.z() * (1.0 - p.z()) * cos(M_PI * p.x());
   //    g(1,3) = -1.75 * p.y() * (1-p.y()) * (1.0 - 2.0 * p.z()) * cos(M_PI * p.x());
   //
   //    g(2,1) = -1.75 *(1-2 * p.x()) * p.z() * (1.0 - p.z()) * cos(M_PI * p.y());
   //    g(2,2) = -1.75 * M_PI * p.x() * (1-p.x()) * p.z() * (1.0 - p.z()) * sin(M_PI * p.x());
   //    g(2,3) = -1.75 * p.x() * (1-p.x()) * (1.0 - 2.0 * p.z()) * cos(M_PI * p.y());
   //
   //    g(3,1) = -0.24 * M_PI * p.z() * p.z() * (1.0 - cos(2.0 * M_PI * p.y())) * sin(2.0*M_PI * p.x());
   //    g(3,2) = -0.24 * M_PI * p.z() * p.z() * (1.0 - cos(2.0 * M_PI * p.x())) * sin(2.0*M_PI * p.y());
   //    g(3,3) = -0.96 * p.z() * std::pow(sin(M_PI * p.x()),2.0) * std::pow(sin(M_PI * p.y()),2.0) + 0.15;
   //    return g;
   // };

   T alpha = 0.2;
   T beta = 0.2;

   auto load = [elas_param](const point<T,3>& p) -> result_type {
      return result_type{0.0,0.0, 0.0};
   };

   auto solution = [elas_param, alpha, beta](const point<T,3>& p) -> result_type {
      T lambda = elas_param.lambda;
      T gamma = alpha + beta + alpha * beta;
      T fx = (1.0/lambda + alpha) * p.x() + /* f(Y)= */ 0.0;
      T fy = (1.0 - gamma/(1.0 + gamma)) * p.y();
      T fz = (1.0/lambda + beta) * p.x() + /* g(X)= */ 0.0 + /* h(Y)= */ 0.0;

      return result_type{fx,fy,fz};
   };

   auto gradient = [elas_param, alpha, beta](const point<T,3>& p) -> result_grad_type {
      T lambda = elas_param.lambda;
      T gamma = alpha + beta + alpha * beta;
      result_grad_type grad = result_grad_type::Zero();

      grad(0,0) = (1.0/lambda + alpha);
      grad(0,1) = /* f'(Y)= */ 0.0;

      grad(1,1) = (1.0 - gamma/(1.0 + gamma));

      grad(2,0) = /* g'(X)= */ 0.0;
      grad(2,1) = /* h'(Y)= */ 0.0;
      grad(2,2) = (1.0/lambda + beta);

      return grad;
   };


   auto neumann = [elas_param](const point<T,3>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.0;
      T fz = 0.0;

      return result_type{fx,fy,fz};
   };

   std::vector<size_t> boundary_neumann(0);


   hyperelasticity_solver<Mesh, T, 3, Storage,  point<T, 3> > nl(msh, rp.degree, elas_param);
   nl.verbose(false);

   nl.compute_initial_state(boundary_neumann);


   auto solve_info = nl.compute(load, solution, neumann, boundary_neumann, rp.n_time_step);

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
void test_hexahedra_diskpp(const run_params& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-2-2-2.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-4-4-4.hex");
   paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-8-8-8.hex");
   paths.push_back(".../diskpp/meshes/3D_hexa/diskpp/testmesh-16-16-16.hex");
   //paths.push_back("../diskpp/meshes/3D_hexa/diskpp/testmesh-32-32-32.hex");

   std::vector<error_type> error_sumup;

   for(int i = 0; i <= runs; i++){
      auto msh = disk::load_cartesian_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_hexahedra_fvca6(const run_params& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_2x2x2.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_4x4x4.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_8x8x8.msh");
   paths.push_back("../diskpp/meshes/3D_hexa/fvca6/hexa_16x16x16.msh");
   //paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_32x32x32.hex");

   std::vector<error_type> error_sumup;

   for(int i = 0; i <= runs; i++){
      auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param));
   }
   printResults(error_sumup);
}

template< typename T>
void test_tetrahedra_netgen(const run_params& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet1.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet2.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet3.mesh");
   paths.push_back("../diskpp/meshes/3D_tetras/netgen/fvca6_tet4.mesh");

   std::vector<error_type> error_sumup;

   for(int i = 0; i <= runs; i++){
      auto msh = disk::load_netgen_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_polyhedra_fvca6(const run_params& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 3;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_10.msh");
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_20.msh");
   paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_30.msh");
   //paths.push_back("../diskpp/meshes/3D_general/fvca6/dbls_40.msh");

   std::vector<error_type> error_sumup;

   for(int i = 0; i <= runs; i++){
      auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param));
   }
   printResults(error_sumup);
}


template< typename T>
void test_tetrahedra_fvca6(const run_params& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.1.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.2.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.3.msh");
   paths.push_back("../diskpp/meshes/3D_tetras/fvca6/tet.4.msh");

   std::vector<error_type> error_sumup;

   for(int i = 0; i <= runs; i++){
      auto msh = disk::load_fvca6_3d_mesh<T>(paths[i].c_str());
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
   rp.verbose  = true;
   rp.n_time_step = 1;

   ElasticityParameters param = ElasticityParameters();
   param.lambda = 1.0;
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

   std::cout << " Test convergence rates in 3D for: "<< std::endl;
   std::cout << " ** Degree k = " << rp.degree << std::endl;
   std::cout << " ** Stab tau = " << param.tau << std::endl;
   std::cout << " ** mu = " << param.mu << std::endl;
   std::cout << " ** lambda = " << param.lambda << std::endl;
   std::cout << " "<< std::endl;

   tc.tic();
   std::cout << "-Tetrahedras fvca6:" << std::endl;
   test_tetrahedra_fvca6<RealType>(rp, param);
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
   test_hexahedra_diskpp<RealType>(rp, param);
   tc.toc();
   std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
   std::cout << " "<< std::endl;


   tc.tic();
   std::cout << "-Polyhedra:"  << std::endl;
   test_polyhedra_fvca6<RealType>(rp, param);
   tc.toc();
   std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
   std::cout << " "<< std::endl;

}
