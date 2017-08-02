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


void
usage(const char *progname)
{
   printf("Usage: %s <options> <filename>\n\n", progname);
   printf("    -k: face degree (>0)\n");
   printf("    -l: difference beetween cell and face degree (-1 <= l <= 1) \n");
   printf("    -n: number of time step (>0)\n");
   printf("    -m: number of sublevel time step (>0)\n");
   printf("    -o: type of neohookean law (1 <= l <= 4) \n");
   printf("    -v: verbose\n");
   printf("    -t: stabilization parameter (>0)\n");
   printf("    -s: adaptative stabilisation\n");
}


template<template<typename, size_t , typename> class Mesh,
         typename T, typename Storage>
void
run_hyperelasticity_solver(const Mesh<T, 2, Storage>& msh, const ParamRun<T>& rp, const ElasticityParameters& elas_param, T gamma)
{
   typedef Mesh<T, 2, Storage> mesh_type;
   typedef static_vector<T, 2> result_type;
   typedef static_matrix<T, 2, 2> result_grad_type;

   std::cout << "Degree k= " << rp.m_face_degree << std::endl;
   std::cout << "mu= " << elas_param.mu << std::endl;
   std::cout << "lambda= " << elas_param.lambda << std::endl;
   std::cout << "gamma= " << gamma << std::endl;
   std::cout << "gamma_normalise= " << gamma/elas_param.mu << std::endl;
   std::cout << "law= " << elas_param.type_law << std::endl;

   auto load = [gamma](const point<T,2>& p) -> result_type {
      T fx = 0.0;
      T fy = gamma;
      return result_type{fx,fy};
   };

   auto solution = [](const point<T,2>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.0;

      return result_type{fx,fy};
   };

   auto gradient = [](const point<T,2>& p) -> result_grad_type {
      result_grad_type grad = result_grad_type::Zero();
      return grad;
   };

   auto neumann = [](const point<T,2>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.0;

      return result_type{fx,fy};
   };

   BoundaryType N1;
   N1.id = 4;//4
   N1.boundary_type = NEUMANN;

   std::vector<BoundaryType> boundary_neumann = {N1}; //by default 0 is for a dirichlet face
   // 4 for Aurrichio test1
   std::vector<BoundaryType> boundary_dirichlet = {};

   hyperelasticity_solver<Mesh, T, 2, Storage,  point<T, 2> >
      nl(msh, rp, elas_param, boundary_neumann, boundary_dirichlet);


   nl.compute_initial_state();


   if(nl.verbose()){
      std::cout << "Solving the problem ..."  << '\n';
   }

   SolverInfo solve_info = nl.compute(load, solution, neumann);

   if(nl.verbose()){
      std::cout << "Total time to solve the problem: " << solve_info.m_time_solver << " sec" << '\n';
   }


   if(nl.test_convergence()){
      std::cout << "l2 error displacement: " << nl.compute_l2_error(solution) << std::endl;
      std::cout << "l2 error gradient: " << nl.compute_l2_gradient_error(gradient) << std::endl;
      std::cout << "Post-processing: " << std::endl;
      nl.compute_discontinuous_solution("sol_disc2D.msh");
      nl.plot_J_at_gausspoint("J_gp2D.msh");
   }
}


int main(int argc, char **argv)
{
    using RealType = double;

    char    *mesh_filename  = nullptr;
    char    *plot_filename  = nullptr;
    int     degree          = 1;

    ParamRun<RealType> rp;

    ElasticityParameters param = ElasticityParameters();

    param.mu = 40.0;
    param.lambda = param.mu * 1E5;
    param.type_law = 3;

    RealType gamma = 1.0;

    int ch;

    while ( (ch = getopt(argc, argv, "g:k:r:v")) != -1 )
   {
      switch(ch)
      {
           case 'g': gamma = atof(optarg); break;
           case 'k':
               degree = atoi(optarg);
               if (degree < 0)
               {
                   std::cout << "Degree must be positive. Falling back to 1." << std::endl;
                   degree = 1;
               }
               rp.m_cell_degree = degree;
               rp.m_face_degree = degree;
               rp.m_grad_degree = degree;
               break;

           case 'r':
                  if(!rp.readParameters(optarg))
                  exit(1);

                  break;

           case 'v':
               rp.m_verbose = true;
               break;

           case '?':
           default:
               std::cout << "wrong arguments" << std::endl;
               exit(1);
      }
   }

    argc -= optind;
    argv += optind;

    if (argc == 0)
    {
        std::cout << "Mesh format: 1D uniform (Not avaible)" << std::endl;
        return 0;
    }

    mesh_filename = argv[0];


    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
        run_hyperelasticity_solver(msh, rp, param, gamma);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad2$") ))
    {
       std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
       auto msh = disk::load_cartesian_2d_mesh2<RealType>(mesh_filename);
       run_hyperelasticity_solver(msh, rp, param, gamma);
       return 0;
    }
}
