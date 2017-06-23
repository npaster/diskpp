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
#include "pLaplaceParameters.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "pLaplace_solver.hpp"

struct run_params
{
    size_t  degree;
    int     l;
    bool    verbose;
};


template<template<typename, size_t , typename> class Mesh,
         typename T, typename Storage>
void
run_pLaplace_solver(const Mesh<T, 2, Storage>& msh, run_params& rp, pLaplaceParameters lap_param)
{
   typedef T result_type;


   auto load = [lap_param](const point<T,2>& pt) -> result_type {
      const size_t p = lap_param.p;
      result_type norm_G = M_PI * sqrt( std::pow(cos(pt.x() * M_PI) * sin(pt.y() * M_PI),2) +
         std::pow(sin(pt.x() * M_PI) * cos(pt.y() * M_PI),2));

      T fx = sin(M_PI*pt.x())*sin(M_PI*pt.y())*
      (
         2.*std::pow(M_PI,2)*std::pow(norm_G, p-2.0) - (p-2)*std::pow(M_PI,4)*std::pow(norm_G, p-4.0)*(
            std::pow(cos(M_PI*pt.x()),2)*cos(2.*M_PI*pt.y()) + cos(2.*M_PI*pt.x())*std::pow(cos(M_PI*pt.y()),2))
      );

      return result_type{fx};
   };

   auto solution = [](const point<T,2>& pt) -> result_type {
      return sin(pt.x() * M_PI) * sin(pt.y() * M_PI);
   };


   pLaplace_solver<Mesh, T, 2, Storage,  point<T, 2> > nl(msh, rp.degree, lap_param);
   nl.verbose(rp.verbose);


   nl.compute_initial_state();

   if(nl.verbose()){
      std::cout << "Solve the problem for p = " << lap_param.p << '\n';
   }

   const size_t n_time_step = 1;

   solve_info solve_info = nl.compute(load, solution, n_time_step);

   if(nl.verbose()){
      std::cout << "Total time to solve the problem: " << solve_info.time_solver << " sec" << '\n';
   }


   if(nl.test_convergence()){

        std::cout << "l2 error: " << nl.compute_l2_error(solution) << std::endl;

        std::cout << "Post-processing: " << std::endl;
        nl.plot_solution_at_gausspoint("sol_2d.msh");
        nl.plot_l2error_at_gausspoint("error_gp_2d.msh", solution);
   }
}

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
void
run_pLaplace_solver(const Mesh<T, 3, Storage>& msh, run_params& rp, const pLaplaceParameters lap_param)
{
   typedef T result_type;

    auto load = [lap_param](const point<T,3>& pt) -> auto {
       T fx = -(std::pow(3.0, lap_param.p/2.0) *(lap_param.p-1)*exp((lap_param.p-1)*(pt.x() + pt.y() + pt.z())));
       return result_type{fx};
   };

   auto solution = [](const point<T,3>& pt) -> auto {
      T fx = exp(pt.x() + pt.y() + pt.z());

      return result_type{fx};
   };


   pLaplace_solver<Mesh, T, 3, Storage,  point<T, 3> > nl(msh, rp.degree, lap_param);
   nl.verbose(rp.verbose);

   nl.compute_initial_state();

   if(nl.verbose()){
      std::cout << "Solve the problem for p = " << lap_param.p << '\n';
   }

   const size_t n_time_step = 1;

   auto solve_info = nl.compute(load, solution, n_time_step);

   if(nl.verbose()){
      std::cout << "Total time to solve the problem: " << solve_info.time_solver << " sec" << '\n';
   }

   if(nl.test_convergence()){

        std::cout << "l2 error: " << nl.compute_l2_error(solution) << std::endl;

        std::cout << "Post-processing: " << std::endl;
        nl.plot_solution_at_gausspoint("sol_3d.msh");
        nl.plot_l2error_at_gausspoint("error_gp_3d.msh", solution);
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

    pLaplaceParameters param = pLaplaceParameters();
    param.p = 2;
    param.tau = 1.0;

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
                elems_1d = atoi(optarg);
                if (elems_1d < 0)
                {
                    std::cout << "Num of elems must be positive. Falling back to 8." << std::endl;
                    elems_1d = 8;
                }
                break;

            case 'p':
                param.p = atoi(optarg);
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

    if (argc == 0)
    {
        std::cout << "Mesh format: 1D uniform (Not avaible)" << std::endl;
        return 0;
    }

    mesh_filename = argv[0];

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<RealType>(mesh_filename);
        run_pLaplace_solver(msh, rp, param);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
        run_pLaplace_solver(msh, rp, param);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<RealType>(mesh_filename);
        run_pLaplace_solver(msh, rp, param);
        return 0;
    }

    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<RealType>(mesh_filename);
        run_pLaplace_solver(msh, rp, param);
        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<RealType>(mesh_filename);
        run_pLaplace_solver(msh, rp, param);
        return 0;
    }

}