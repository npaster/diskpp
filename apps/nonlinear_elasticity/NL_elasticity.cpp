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

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "NL_elasticity_solver.hpp"

struct run_params
{
    size_t  degree;
    int     l;
    bool    verbose;
};

// template<template<typename, size_t , typename> class Mesh,
// typename T, typename Storage>
// void
// run_NL_elasticity_solver(const Mesh<T, 1, Storage>& msh, run_params& rp)
// {
//    typedef Mesh<T, 1, Storage> mesh_type;
//
// //    auto load = [](const point<T,1>& p) -> auto {
// //       const T lambda = T{0.0};
// //       const T mu = T{1.0};
// //       return (T{2.0} * mu + lambda) * M_PI *  M_PI * sin(p.x() * M_PI);
// //    };
// //
// //    auto solution = [](const point<T,1>& p) -> auto {
// //       return sin(p.x() * M_PI);
// //    };
//
//    auto load = [](const point<T,1>& p) -> auto {
//       const T lambda = T{0.0}3
//       const T mu = T{1.0};
//       return 0.0;
//    };
//
//    auto solution = [](const point<T,1>& p) -> auto {
//       return p.x();
//    };
//
//    NL_elasticity_solver<mesh_type,  point<T, 1> > nl(msh, rp.degree);
//    nl.verbose(rp.verbose);
//
//    //auto info_offline = nl.compute_offline();
//
//    //    if(nl.verbose()){
//    //       std::cout << "Off_line computations: " << info_offline.time_offline << " sec"  << '\n';
//    //    }
//
//    nl.compute_initial_state();
//
//    if(nl.verbose()){
//       std::cout << "Solve the problem: " << '\n';
//    }
//
//    const size_t n_time_step = 1;
//
//    solve_info solve_info = nl.compute(load, solution, n_time_step);
//
//    if(nl.verbose()){
//       std::cout << "Total time to solve the problem: " << solve_info.time_solver << " sec" << '\n';
//    }
//
//
//    if(nl.test_convergence()){
//
//       std::cout << "l2 error: " << nl.compute_l2_error(solution) << std::endl;
//
//       std::cout << "Post-processing: " << std::endl;
//       //nl.plot_solution_at_gausspoint("sol_elas_2d.msh");
//       //nl.plot_l2error_at_gausspoint("error_gp_2d_.msh", solution);
//    }
// }

template<template<typename, size_t , typename> class Mesh,
         typename T, typename Storage>
void
run_NL_elasticity_solver(const Mesh<T, 2, Storage>& msh, run_params& rp)
{
   typedef Mesh<T, 2, Storage> mesh_type;
   typedef static_vector<T, 2> result_type;

   auto loadg = [](const point<T,2>& p) -> result_type {
      const T lambda =0.0;
      T fx = 4.*M_PI*M_PI*sin(M_PI*p.x())*sin(M_PI*p.y());
      T fy = 4.*M_PI*M_PI*cos(M_PI*p.x())*cos(M_PI*p.y());
      return result_type{fx,fy};
   };

   auto solutiong = [](const point<T,2>& p) -> result_type {
      const T lambda =1.0;
      T fx = sin(M_PI*p.x())*sin(M_PI*p.y());
      T fy = cos(M_PI*p.x())*cos(M_PI*p.y());

      return result_type{fx,fy};
   };


   auto load_square = [](const point<T,2>& p) -> result_type {
      const T lambda =0.0;
      const T mu = 1.0;
      T scale = 0.01;
      T fx = scale * scale * 4.*(lambda + 2.0 * mu) * (2.0 *p.x() +1.0);
      T fy = scale * scale * 4.*(lambda + 2.0 * mu) * (2.0 *p.y() +1.0);
      return result_type{fx,fy};
   };

   auto solution_square = [](const point<T,2>& p) -> result_type {
      T scale = 0.01;
      T fx = scale * p.x() * p.x();
      T fy = scale * p.y() * p.y();

      return result_type{fx,fy};
   };



   auto load_lin = [](const point<T,2>& p) -> result_type {
      return result_type{0.0,0.0};
   };

   auto solution_lin = [](const point<T,2>& p) -> result_type {
      const T scale = -0.1;
      T fx = scale * p.x();
      T fy = scale * p.y();

      return result_type{fx,fy};
   };

   auto load_traction = [](const point<T,2>& p) -> result_type {
      return result_type{0.0,0.0};
   };

   auto solution_traction = [](const point<T,2>& p) -> result_type {
      const T scale = 0.1;
      T fx = 0.0;
      T fy = scale * p.y();

      return result_type{fx,fy};
   };

   NL_elasticity_solver<Mesh, T, 2, Storage,  point<T, 2> > nl(msh, rp.degree);
   nl.verbose(rp.verbose);

   auto info_offline = nl.compute_offline();

   if(nl.verbose()){
      std::cout << "Off_line computations: " << info_offline.time_offline << " sec"  << '\n';
   }

   nl.compute_initial_state();

   if(nl.verbose()){
      std::cout << "Solve the problem: " << '\n';
   }

   const size_t n_time_step = 0;

   solve_info solve_info = nl.compute(load_lin, solution_lin, n_time_step);

   if(nl.verbose()){
      std::cout << "Total time to solve the problem: " << solve_info.time_solver << " sec" << '\n';
   }


   if(nl.test_convergence()){

        std::cout << "l2 error: " << nl.compute_l2_error(solution_lin) << std::endl;

        std::cout << "Post-processing: " << std::endl;
        nl.plot_solution_at_gausspoint("sol_elas_2d.msh");
        nl.plot_l2error_at_gausspoint("error_gp_2d_.msh", solution_lin);
   }

   nl.saveMesh2D("msh2d.msh");
   nl.compute_deformed("deforme2d.msh");
   nl.compute_discontinuous_solution("depl2d.msh");
}

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
void
run_NL_elasticity_solver(const Mesh<T, 3, Storage>& msh, run_params& rp)
{
   typedef Mesh<T, 3, Storage> mesh_type;
   typedef static_vector<T, 3> result_type;

    auto load = [](const point<T,3>& p) -> auto {
      const T lambda =1.0;
      T fx = M_PI*M_PI*(12*lambda*cos(2*M_PI*p.x())*sin(2*M_PI*p.y())*sin(2*M_PI*p.z())
      + 8 * (3*cos(2*M_PI*p.x())-1)*sin(2*M_PI*p.y())*sin(2*M_PI*p.z())
      - cos(M_PI*p.x())*sin(M_PI*(p.y()+p.z()))
      + (1 + 3./(1+lambda))* sin(M_PI*p.x())*sin(M_PI*p.y())*sin(M_PI*p.z()) );

      T fy = M_PI*M_PI*(12*lambda*cos(2*M_PI*p.y())*sin(2*M_PI*p.x())*sin(2*M_PI*p.z())
      + 8 * (3*cos(2*M_PI*p.y())-1)*sin(2*M_PI*p.x())*sin(2*M_PI*p.z())
      - cos(M_PI*p.y())*sin(M_PI*(p.x()+p.z())) +
      (1 + 3./(1+lambda))* sin(M_PI*p.x())*sin(M_PI*p.y())*sin(M_PI*p.z()) );

      T fz = M_PI*M_PI*(12*lambda*cos(2*M_PI*p.z())*sin(2*M_PI*p.y())*sin(2*M_PI*p.x())
      + 8 * (3*cos(2*M_PI*p.z())-1)*sin(2*M_PI*p.y())*sin(2*M_PI*p.x())
      - cos(M_PI*p.z())*sin(M_PI*(p.y()+p.x())) +
      (1 + 3./(1+lambda))* sin(M_PI*p.x())*sin(M_PI*p.y())*sin(M_PI*p.z()) );

      return result_type{fx,fy,fz};
   };

   auto solution = [](const point<T,3>& p) -> auto {
      const T lambda =1.0;
      T fx = sin(2*M_PI*p.y())*sin(2*M_PI*p.z())*(-1 + cos(2*M_PI*p.x()))
      + (1./(1+lambda))* sin(M_PI*p.x())*sin(M_PI*p.y())*sin(M_PI*p.z());
      T fy = sin(2*M_PI*p.z())*sin(2*M_PI*p.x())*(-1 + cos(2*M_PI*p.y()))
      + (1./(1+lambda))* sin(M_PI*p.x())*sin(M_PI*p.y())*sin(M_PI*p.z());
      T fz = sin(2*M_PI*p.x())*sin(2*M_PI*p.y())*(-1 + cos(2*M_PI*p.z()))
      + (1./(1+lambda))* sin(M_PI*p.x())*sin(M_PI*p.y())*sin(M_PI*p.z());
      return result_type{fx,fy,fz};
   };


   auto load_lin = [](const point<T,3>& p) -> result_type {
      return result_type{0.0,0.0,0.0};
   };

   auto solution_lin = [](const point<T,3>& p) -> result_type {
      const T scale = 0.1;
      T fx = scale * p.x();
      T fy = scale * p.y();
      T fz = scale * p.z();
      return result_type{fx,fy,fz};
   };


   NL_elasticity_solver<Mesh, T, 3, Storage,  point<T, 3> > nl(msh, rp.degree);
   nl.verbose(rp.verbose);

   auto info_offline = nl.compute_offline();

   if(nl.verbose()){
      std::cout << "Off_line computations: " << info_offline.time_offline << " sec"  << '\n';
   }

   nl.compute_initial_state();

   if(nl.verbose()){
      std::cout << "Solve the problem: " << '\n';
   }

   const size_t n_time_step = 0;

   auto solve_info = nl.compute(load_lin, solution_lin, n_time_step);

   if(nl.verbose()){
      std::cout << "Total time to solve the problem: " << solve_info.time_solver << " sec" << '\n';
   }

   if(nl.test_convergence()){

        std::cout << "l2 error: " << nl.compute_l2_error(solution_lin) << std::endl;

        std::cout << "Post-processing: " << std::endl;
        nl.plot_solution_at_gausspoint("sol_elas_3d.msh");
        nl.plot_l2error_at_gausspoint("error_gp_3d_.msh", solution_lin);
   }

   nl.saveMesh3D("msh3d.msh");
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

    int ch;

    while ( (ch = getopt(argc, argv, "k:l:n:p:v")) != -1 )
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
                plot_filename = optarg;
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
        auto msh = disk::load_uniform_1d_mesh<RealType>(0, 1, elems_1d);
        //run_NL_elasticity_solver(msh, rp);
        return 0;
    }

    mesh_filename = argv[0];

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        auto msh = disk::load_fvca5_2d_mesh<RealType>(mesh_filename);
        run_NL_elasticity_solver(msh, rp);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
        run_NL_elasticity_solver(msh, rp);
        return 0;
    }

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<RealType>(mesh_filename);
        run_NL_elasticity_solver(msh, rp);
        return 0;
    }

    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        auto msh = disk::load_netgen_3d_mesh<RealType>(mesh_filename);
        run_NL_elasticity_solver(msh, rp);
        return 0;
    }

    /* DiSk++ cartesian 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
        auto msh = disk::load_cartesian_3d_mesh<RealType>(mesh_filename);
        run_NL_elasticity_solver(msh, rp);
        return 0;
    }

}
