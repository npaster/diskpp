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

#include "hyperelasticity2_solver.hpp"


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

   std::cout << "Degree k= " << rp.m_degree << std::endl;
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

   BoundaryConditions N1;
   N1.id = 4;//4
   N1.boundary_type = NEUMANN;

   std::vector<BoundaryConditions> boundary_neumann = {N1}; //by default 0 is for a dirichlet face
   // 4 for Aurrichio test1
   std::vector<BoundaryConditions> boundary_dirichlet = {};

   hyperelasticity2_solver<Mesh, T, 2, Storage,  point<T, 2> > nl(msh, rp, elas_param);


   nl.compute_initial_state(boundary_neumann, boundary_dirichlet);


   if(nl.verbose()){
      std::cout << "Solving the problem ..."  << '\n';
   }

   SolverInfo solve_info = nl.compute(load, solution, neumann, boundary_neumann, boundary_dirichlet);

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
    int     degree          = 1;
    int     l               = 0;
    int     n_time_step     = 1;
    int     sublevel        = 1;

    ParamRun<RealType> rp;
    rp.m_sublevel = 4;
    rp.m_verbose = false;

    ElasticityParameters param = ElasticityParameters();

    param.mu = 40.0;
    param.lambda = param.mu * 1E5;
    param.type_law = 3;

    RealType gamma = 1.0;

    int ch;

    while ( (ch = getopt(argc, argv, "e:i:g:k:l:m:n:t:v")) != -1 )
    {
       switch(ch)
       {
          case 'i': rp.m_iter_max = atof(optarg); break;
          case 'e': rp.m_epsilon = atof(optarg); break;
          case 'g': gamma = atof(optarg); break;
          case 'k':
             degree = atoi(optarg);
             if (degree < 1)
             {
                std::cout << "Degree must be positive. Falling back to 1." << std::endl;
                degree = 1;
             }
             rp.m_degree = degree;
             break;

          case 'l':
             l = atoi(optarg);
             if (l < -1 or l > 1)
             {
                std::cout << "l can be -1, 0 or 1. Falling back to 0." << std::endl;
                l = 0;
             }
             rp.m_l = l;
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

          case 'n':
             n_time_step = atoi(optarg);
             if (n_time_step <= 0)
             {
                std::cout << "Number of time step must be positive. Falling back to 1." << std::endl;
                n_time_step = 1;
             }
             rp.m_n_time_step = n_time_step;
             break;

          case 't': param.type_law = atoi(optarg); break;

          case 'v': rp.m_verbose = true; break;

          case 'h':
          case '?':
          default:
             std::cout << "wrong arguments" << std::endl;
             usage(argv[0]);
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
