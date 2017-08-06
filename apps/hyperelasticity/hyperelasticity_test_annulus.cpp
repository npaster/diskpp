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
#include "BoundaryConditions.hpp"
#include "Parameters.hpp"
#include "Informations.hpp"

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
run_hyperelasticity_solver(const Mesh<T, 2, Storage>& msh, const ParamRun<T>& rp, const ElasticityParameters& elas_param,
                           const std::string& file_error, const size_t& n=0)
{
   typedef Mesh<T, 2, Storage> mesh_type;
   typedef static_vector<T, 2> result_type;
   typedef static_matrix<T, 2, 2> result_grad_type;

   T r0 = 1.5; T R0 = 0.5;
   T alpha = (r0 - R0)/R0;

   auto load = [](const point<T,2>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.0;
      return result_type{fx,fy};
   };

   auto solution = [alpha](const point<T,2>& p) -> result_type {
      T fx = alpha * p.x();
      T fy = alpha * p.y();

      return result_type{fx,fy};
   };

   auto neumann = [](const point<T,2>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.0;

      return result_type{fx,fy};
   };

   BoundaryType N1;
   N1.id = 4;
   N1.boundary_type = FREE;

   std::vector<BoundaryType> boundary_neumann = {N1};
   std::vector<BoundaryType> boundary_dirichlet = {};

   hyperelasticity_solver<Mesh, T, 2, Storage,  point<T, 2> >
      nl(msh, rp, elas_param, boundary_neumann, boundary_dirichlet);

   nl.compute_initial_state();

   SolverInfo solve_info = nl.compute(load, solution, neumann);

   error_type error;
   error.h = average_diameter(msh);
   error.degree = rp.m_face_degree;
   error.nb_dof = nl.getDofs();
   error.error_depl = 10E6;
   error.error_grad = 10E6;

   if(nl.test_convergence()){
      auto error_annulus = nl.compute_l2_error_annulus(file_error);
      error.error_depl = error_annulus.first;
      error.error_grad = error_annulus.second;
      std::cout << "Post-processing: " << std::endl;
      std::string name = "Result_annulus_k" + std::to_string(rp.m_cell_degree) + "_l" + std::to_string(rp.m_grad_degree)
                        + "_m" + std::to_string(n) + "_" + std::to_string(elas_param.lambda) + "_";
      nl.compute_discontinuous_solution(name + "sol_disc2D.msh");
      nl.compute_conforme_solution(name +"sol_conf2D.msh");
      nl.compute_deformed(name +"def2D.msh");
      nl.plot_displacement_at_gausspoint(name +"depl_gp2D.msh");
      nl.plot_J_at_gausspoint(name +"J_gp2D.msh");
      nl.plot_J(name +"J_dis2d.msh");
      nl.compute_discontinuous_PK1(name +"PK1_disc2D.msh");
      nl.compute_discontinuous_Prr(name +"Prr.msh", "Prr");
      nl.compute_discontinuous_Prr(name +"Poo.msh", "Poo");
      nl.compute_discontinuous_VMIS(name +"VM.msh");
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
void test_annulus(const ParamRun<T>& rp, const ElasticityParameters& elas_param, const std::string& file_error)
{
   size_t runs = 4;

   std::vector<std::string> paths;
   paths.push_back("../meshes/Anneau/anneau_1.medit2d");
   paths.push_back("../meshes/Anneau/anneau_2.medit2d");
   paths.push_back("../meshes/Anneau/anneau_3.medit2d");
   paths.push_back("../meshes/Anneau/anneau_4.medit2d");
   paths.push_back("../meshes/Anneau/anneau_5.medit2d");

   std::vector<error_type> error_sumup;

   for(size_t i = 0; i < runs; i++){
      auto msh = disk::load_medit_2d_mesh<T>(paths[i].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param, file_error, i+1));
   }
   printResults(error_sumup);
}


template< typename T>
void test_annulus_locking(const ParamRun<T>& rp,  ElasticityParameters& elas_param, const std::string& file_error)
{
   size_t runs = 6;

   std::vector<std::string> paths;
   paths.push_back("../meshes/Anneau/anneau_3.medit2d");

   std::array<T,5> lambda_test = {1.66644, 16.6644, 166.644, 1666.44, 16664.4, 166644.0};

   std::vector<error_type> error_sumup;

   for(size_t i = 0; i < runs; i++){
      elas_param.lambda = lambda_test[i];
      auto msh = disk::load_medit_2d_mesh<T>(paths[0].c_str());
      error_sumup.push_back(run_hyperelasticity_solver(msh, rp, elas_param, file_error, 3));
   }
   printResults(error_sumup);
}




int main(int argc, char **argv)
{
   using RealType = double;

   char    *mesh_filename  = nullptr;
   char    *plot_filename  = nullptr;
   std::string    file_error;
   size_t num = 0;
   bool    convergence_rates = false;

   ParamRun<RealType> rp;

   ElasticityParameters param = ElasticityParameters();
   param.lambda = 16664.4;
   param.mu = 0.333;
   param.type_law = 1;

   int ch;

   while ( (ch = getopt(argc, argv, "ce:l:n:r:")) != -1 )
   {
      switch(ch)
      {
         case 'c': convergence_rates = true; break;

         case 'l': param.lambda = atof(optarg); break;

         case 'e': file_error = optarg; break;

         case 'n': num = atoi(optarg); break;

         case 'r':
            if(!rp.readParameters(optarg))
            exit(1);

            break;

         case '?':
         default:
            std::cout << "wrong arguments" << std::endl;
            usage(argv[0]);
            exit(1);
      }
   }

   argc -= optind;
   argv += optind;

    mesh_filename = argv[0];

   timecounter tc;

   std::cout << " Test convergence rates for: "<< std::endl;
   std::cout << " ** Face_Degree = " << rp.m_face_degree << std::endl;
   std::cout << " ** Cell_Degree  = " << rp.m_cell_degree << std::endl;
   std::cout << " ** Grad_Degree  = " << rp.m_grad_degree << std::endl;
   std::cout << " ** mu = " << param.mu << std::endl;
   std::cout << " ** lambda = " << param.lambda << std::endl;
   std::cout << " "<< std::endl;

   if(convergence_rates){
      tc.tic();
      std::cout << "-Medit 2D:" << std::endl;
      test_annulus<RealType>(rp, param, file_error);
      tc.toc();
      std::cout << "Time to test convergence rates: " << tc.to_double() << std::endl;
      std::cout << " "<< std::endl;
   }
   else{
      /* Medit 2d*/
      if (std::regex_match(mesh_filename, std::regex(".*\\.medit2d$") ))
      {
         std::cout << "Guessed mesh format: Medit format" << std::endl;
         auto msh = disk::load_medit_2d_mesh<RealType>(mesh_filename);
         run_hyperelasticity_solver(msh, rp, param, file_error, num);
         return 0;
      }
   }

}
