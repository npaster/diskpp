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
};


template<template<typename, size_t , typename> class Mesh,
         typename T, typename Storage>
void
run_hyperelasticity_solver(const Mesh<T, 2, Storage>& msh, run_params& rp, ElasticityParameters elas_param)
{
   typedef Mesh<T, 2, Storage> mesh_type;
   typedef static_vector<T, 2> result_type;


//    auto load = [](const point<T,2>& p) -> result_type {
//        T fx = -16.0 * p.x() * p.y() * p.y() -4.0 * p.x() * p.x() * p.x();
//        T fy = -16.0 * p.x() * p.x() * p.y() -4.0 * p.y() * p.y() * p.y();
//       return result_type{fx,fy};
//    };
//
//    auto solution = [](const point<T,2>& p) -> result_type {
//       T fx = 1.2 * p.x() * p.x() * p.y();
//       T fy = (p.x() * p.y() * p.y()) * sin(p.x()) ;
//
//       return result_type{fx,fy};
//    };


   auto load = [elas_param](const point<T,2>& p) -> result_type {
      T lambda = elas_param.lambda;
      T mu = elas_param.mu;

      T num11 = lambda*(18 * p.x() * std::pow(p.y(), 2.0) + 6.0 * p.x()) + 6.0 * p.x();
      T num21 = lambda*(18 * p.y() * std::pow(p.x(), 2.0) + 6.0 * p.y()) + 6.0 * p.y();

      T dem11 = 9.0 * std::pow(p.x(), 2.0) * std::pow(p.y(), 2.0) + 3.0 * (std::pow(p.x(), 2.0) + std::pow(p.y(), 2.0)) * ( 1.0 + 1.0/lambda)
      + 1.0 + 2.0 / lambda + 1.0/std::pow(lambda, 2.0);

      T dem12 = 3.0 * std::pow(p.x(), 2.0) + 1.0 + 1.0 / lambda;
      T dem21 = 3.0 * std::pow(p.y(), 2.0) + 1.0 + 1.0 / lambda;

      T fx = - num11 / (dem12 *dem11) - 6.0*mu*p.x() + 6.0*p.x()*(lambda * log(dem11) - mu)/std::pow(dem12,2.0);
      T fy = - num21 / (dem21 *dem11) - 6.0*mu*p.y() + 6.0*p.y()*(lambda * log(dem11) - mu)/std::pow(dem21,2.0);
      return result_type{fx,fy};
   };

   auto solution = [elas_param](const point<T,2>& p) -> result_type {
      T lambda = elas_param.lambda;
      T fx = std::pow(p.x(), 3.0) + p.x()/lambda;
      T fy = std::pow(p.y(), 3.0) + p.y()/lambda;

      return result_type{fx,fy};
   };


   hyperelasticity_solver<Mesh, T, 2, Storage,  point<T, 2> > nl(msh, rp.degree, elas_param);
   nl.verbose(rp.verbose);


   nl.compute_initial_state();

   if(nl.verbose()){
      std::cout << "Solving the problem ..."  << '\n';
   }

   solve_info solve_info = nl.compute(load, solution, rp.n_time_step);

   if(nl.verbose()){
      std::cout << "Total time to solve the problem: " << solve_info.time_solver << " sec" << '\n';
   }


   if(nl.test_convergence()){
      std::cout << "avrage diameter h: " << average_diameter(msh) << std::endl;
      std::cout << "l2 error: " << nl.compute_l2_error(solution) << std::endl;

        std::cout << "Post-processing: " << std::endl;
        nl.compute_discontinuous_solution("sol_disc2D.msh");
        nl.compute_deformed("def2D.msh");
   }
}

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
void
run_hyperelasticity_solver(const Mesh<T, 3, Storage>& msh, run_params& rp, const ElasticityParameters elas_param)
{
   typedef Mesh<T, 3, Storage> mesh_type;
   typedef static_vector<T, 3> result_type;

   auto load = [](const point<T,3>& p) -> result_type {
      return result_type{0.0,0.0,0.0};
   };

   auto solution = [](const point<T,3>& p) -> result_type {
      T fx = -1.75 * p.y() * (1-p.y()) * p.z() * (1.0 - p.z()) * cos(M_PI * p.x());
      T fy = -1.75 * p.x() * (1-p.x()) * p.z() * (1.0 - p.z()) * cos(M_PI * p.y());
      T fz = -0.12 * p.z() * p.z() * (1.0 - cos(2.0 * M_PI * p.x())) * (1.0 - cos(2.0*M_PI * p.y())) + 0.15 * p.z();
      return result_type{fx,fy,fz};
   };


   auto gradient = [](const point<T,3>& p) -> static_matrix<T, 3, 3> {
      static_matrix<T, 3, 3> g;
      g(1,1) = -1.75 * M_PI * p.y() * (1-p.y()) * p.z() * (1.0 - p.z()) * sin(M_PI * p.x());
      g(1,2) = -1.75 *(1-2 * p.y()) * p.z() * (1.0 - p.z()) * cos(M_PI * p.x());
      g(1,3) = -1.75 * p.y() * (1-p.y()) * (1.0 - 2.0 * p.z()) * cos(M_PI * p.x());

      g(2,1) = -1.75 *(1-2 * p.x()) * p.z() * (1.0 - p.z()) * cos(M_PI * p.y());
      g(2,2) = -1.75 * M_PI * p.x() * (1-p.x()) * p.z() * (1.0 - p.z()) * sin(M_PI * p.x());
      g(2,3) = -1.75 * p.x() * (1-p.x()) * (1.0 - 2.0 * p.z()) * cos(M_PI * p.y());

      g(3,1) = -0.24 * M_PI * p.z() * p.z() * (1.0 - cos(2.0 * M_PI * p.y())) * sin(2.0*M_PI * p.x());
      g(3,2) = -0.24 * M_PI * p.z() * p.z() * (1.0 - cos(2.0 * M_PI * p.x())) * sin(2.0*M_PI * p.y());
      g(3,3) = -0.96 * p.z() * std::pow(sin(M_PI * p.x()),2.0) * std::pow(sin(M_PI * p.y()),2.0) + 0.15;
      return g;
   };


   hyperelasticity_solver<Mesh, T, 3, Storage,  point<T, 3> > nl(msh, rp.degree, elas_param);
   nl.verbose(rp.verbose);

   nl.compute_initial_state();

   if(nl.verbose()){
      std::cout << "Solving the problem ..." << '\n';
   }

   auto solve_info = nl.compute(load, solution, rp.n_time_step);

   if(nl.verbose()){
      std::cout << "Total time to solve the problem: " << solve_info.time_solver << " sec" << '\n';
   }

   if(nl.test_convergence()){
      std::cout << "avrage diameter h: " << average_diameter(msh) << std::endl;
      std::cout << "l2 error: " << nl.compute_l2_error(solution) << std::endl;

        std::cout << "Post-processing: " << std::endl;
        nl.compute_discontinuous_solution("sol_disc3D.msh");
        nl.compute_deformed("def3D.msh");
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

    ElasticityParameters param = ElasticityParameters();
    param.lambda = 1.0;
    param.mu = 1.0;
    param.tau = 1000.0;
    param.adaptative_stab = false;

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
        run_hyperelasticity_solver(msh, rp, param);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
        run_hyperelasticity_solver(msh, rp, param);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<RealType>(mesh_filename);
        run_hyperelasticity_solver(msh, rp, param);
        return 0;
    }

    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<RealType>(mesh_filename);
        run_hyperelasticity_solver(msh, rp, param);
        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<RealType>(mesh_filename);
        run_hyperelasticity_solver(msh, rp, param);
        return 0;
    }

}
