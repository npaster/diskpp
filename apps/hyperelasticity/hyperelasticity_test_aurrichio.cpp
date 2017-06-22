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

struct run_params
{
    size_t  degree;
    int     l;
    bool    verbose;
    size_t n_time_step;
    double gamma;
};


template<template<typename, size_t , typename> class Mesh,
         typename T, typename Storage>
void
run_hyperelasticity_solver(const Mesh<T, 2, Storage>& msh, const run_params& rp, const ElasticityParameters& elas_param)
{
   typedef Mesh<T, 2, Storage> mesh_type;
   typedef static_vector<T, 2> result_type;
   
   std::cout << "Degree k= " << rp.degree << std::endl;
   std::cout << "mu= " << elas_param.mu << std::endl;
   std::cout << "lambda= " << elas_param.lambda << std::endl;
   std::cout << "tau= " << elas_param.tau << std::endl;
   std::cout << "gamma= " << rp.gamma << std::endl;
   std::cout << "gamma_normalise= " << rp.gamma/elas_param.mu << std::endl;
   std::cout << "law= " << elas_param.type_law << std::endl;
   
   auto load = [rp](const point<T,2>& p) -> result_type {
      T fx = 0.0;
      T fy = rp.gamma;
      return result_type{fx,fy};
   };
   
   auto solution = [elas_param](const point<T,2>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.0;
      
      return result_type{fx,fy};
   };

   auto neumann = [elas_param](const point<T,2>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.0;
      
      return result_type{fx,fy};
   };

   std::vector<size_t> boundary_neumann(1,4); //by default 0 is for a dirichlet face
   // 4 for Aurrichio test1 

   hyperelasticity_solver<Mesh, T, 2, Storage,  point<T, 2> > nl(msh, rp.degree, elas_param);
   nl.verbose(rp.verbose);


   nl.compute_initial_state();

   if(nl.verbose()){
      std::cout << "Solving the problem ..."  << '\n';
   }

   solve_info solve_info = nl.compute(load, solution, neumann, boundary_neumann, rp.n_time_step);

   if(nl.verbose()){
      std::cout << "Total time to solve the problem: " << solve_info.time_solver << " sec" << '\n';
   }


   if(nl.test_convergence()){
        std::cout << "Post-processing: " << std::endl;
        nl.compute_discontinuous_solution("sol_disc2D.msh");
        nl.compute_deformed("def2D.msh");
   }
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
    rp.gamma = 200.0;

    ElasticityParameters param = ElasticityParameters();
    
    param.mu = 40.0;
    param.lambda = param.mu * 10E5;
    param.tau = 10.0;
    param.adaptative_stab = false;
    param.type_law = 1;

    int ch;

    while ( (ch = getopt(argc, argv, "k:l:n:p:v:t:g:j")) != -1 )
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
               
            case 'g':
               rp.gamma = atof(optarg);
               break;
               
            case 'j':
               param.type_law = atoi(optarg);
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
        run_hyperelasticity_solver(msh, rp, param);
        return 0;
    }
    
    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad2$") ))
    {
       std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
       auto msh = disk::load_cartesian_2d_mesh2<RealType>(mesh_filename);
       run_hyperelasticity_solver(msh, rp, param);
       return 0;
    }
}
