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



template<template<typename, size_t , typename> class Mesh,
         typename T, typename Storage>
void
run_hyperelasticity_solver(const Mesh<T, 2, Storage>& msh, ParamRun<T>& rp, const ElasticityParameters& elas_param)
{
   typedef static_vector<T, 2> result_type;
   typedef static_matrix<T, 2, 2> result_grad_type;


T alpha = 0.3;

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


//    T r0 = 0.6; T R0 = 0.5;
//    T alpha = 1.0;// (r0 - R0)/R0;
//    auto load = [elas_param, alpha](const point<T,2>& p) -> result_type {
//       T lambda = elas_param.lambda;
//       T mu = elas_param.mu;
//
//       T fx = 0.0;
//       T fy = 0.0;
//       return result_type{fx,fy};
//    };
//
//    auto solution = [elas_param, alpha](const point<T,2>& p) -> result_type {
//       T lambda = elas_param.lambda;
//       T fx = alpha *  p.x();
//       T fy = alpha *  p.y();
//
//       return result_type{fx,fy};
//    };
//
//    auto gradient = [elas_param, alpha](const point<T,2>& p) -> result_grad_type {
//       T lambda = elas_param.lambda;
//       result_grad_type grad = result_grad_type::Zero();
//       return grad;
//    };



   auto neumann = [elas_param](const point<T,2>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.0;

      return result_type{fx,fy};
   };

//    BoundaryType N1;
//    N1.id = 7;//4,11
//    N1.boundary_type = FREE;
//
//
//    BoundaryType N2;
//    N2.id = 10;//14
//    N2.boundary_type = FREE;
//
//
//    BoundaryType N4;
//    N4.id = 4;//14
//    N4.boundary_type = NEUMANN;

   const std::vector<BoundaryType> boundary_neumann = {}; //by default 0 is for a dirichlet face
   // 4 for Aurrichio test1


//    BoundaryType d3;
//    d3.id = 3;
//    d3.boundary_type = CLAMPED;
//
//
//    BoundaryType d1;
//    d1.id = 1;
//    d1.boundary_type = DX;

   const std::vector<BoundaryType> boundary_dirichlet = {};

   // Solve
   hyperelasticity_solver<Mesh, T, 2, Storage,  point<T, 2> >
      nl(msh, rp, elas_param, boundary_neumann, boundary_dirichlet);


   nl.compute_initial_state();

   if(nl.verbose()){
      std::cout << "Solving the problem ..."  << '\n';
   }

   const SolverInfo solve_info = nl.compute(load, solution, neumann);

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

   // PostProcessing
   if(nl.test_convergence()){
      std::cout << "average diameter h: " << average_diameter(msh) << std::endl;
      std::cout << "l2 error displacement: " << nl.compute_l2_error_displacement(solution) << std::endl;
      std::cout << "l2 error gradient: " << nl.compute_l2_error_gradient(gradient) << std::endl;
      try{
         std::cout << "l2 error energy: " << nl.compute_l2_error_energy(gradient) << std::endl;
      }
      catch(const std::invalid_argument& ia){
         std::cerr << "Invalid argument: " << ia.what() << " in compute_l2_error_energy" << std::endl;
      }
      try{
         std::cout << "l2 error PK1: " << nl.compute_l2_error_PK1(gradient) << std::endl;
      }
      catch(const std::invalid_argument& ia){
         std::cerr << "Invalid argument: " << ia.what() << " in compute_l2_error_PK1" << std::endl;
      }
        std::cout << "Post-processing: " << std::endl;
        nl.compute_discontinuous_displacement("depl_disc2D.msh");
        nl.compute_continuous_displacement("depl_cont2D.msh");
        nl.compute_deformed("def2D.msh");
        nl.compute_J_GP("J_GP.msh");
        nl.compute_continuous_J("J_cont.msh");
        nl.compute_discontinuous_J("J_disc.msh");
        try {
           nl.compute_discontinuous_VMIS("VM_disc.msh");
        }
        catch(const std::invalid_argument& ia){
           std::cerr << "Invalid argument: " << ia.what() << " in VMIS_disc" << std::endl;
        }
        try {
           nl.compute_continuous_VMIS("VM_cont.msh");
        }
        catch(const std::invalid_argument& ia){
           std::cerr << "Invalid argument: " << ia.what() << " in VMIS_disc" << std::endl;
        }
        try {
           nl.compute_VMIS_GP("VM_GP.msh");
        }
        catch(const std::invalid_argument& ia){
           std::cerr << "Invalid argument: " << ia.what() << " in VMIS_disc" << std::endl;
        }
        try {
           nl.compute_l2error_VMIS_GP("VM_GP_error.msh", gradient);
        }
        catch(const std::invalid_argument& ia){
           std::cerr << "Invalid argument: " << ia.what() << " in VMIS_disc" << std::endl;
        }
        nl.plot_analytical_VMIS("VM_true.msh", gradient);
   }
}

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
void
run_hyperelasticity_solver(const Mesh<T, 3, Storage>& msh, ParamRun<T>& rp, const ElasticityParameters elas_param)
{
   typedef static_vector<T, 3> result_type;
   typedef static_matrix<T, 3, 3> result_grad_type;


//    T alpha = 0.2;
//    T beta = 0.2;
//    T factor = 0.5;
//
//
//    auto load = [elas_param, alpha,beta, factor](const point<T,3>& p) -> result_type {
//       T fx = M_PI*M_PI*elas_param.mu*alpha * sin(M_PI*p.y());
//       T fy = 0.0;
//       T fz = M_PI*M_PI*elas_param.mu*beta * sin(M_PI*p.x());
//       return factor*result_type{fx, fy, fz};
//    };
//
//    auto solution = [elas_param, alpha, beta, factor](const point<T,3>& p) -> result_type {
//       T lambda = elas_param.lambda;
//       T gamma = alpha + beta + alpha * beta;
//       T fx = (1.0/lambda + alpha) * p.x() + /* f(Y)= */ alpha * sin(M_PI*p.y());
//       T fy = (1.0/lambda - gamma/(1.0 + gamma)) * p.y();
//       T fz = (1.0/lambda + beta) * p.z() + /* g(X)= */ beta * sin(M_PI*p.x()) + /* h(Y)= */ 0.0;
//
//       return factor*result_type{fx,fy,fz};
//    };
//
//    auto gradient = [elas_param, alpha, beta, factor](const point<T,3>& p) -> result_grad_type {
//       T lambda = elas_param.lambda;
//       T gamma = alpha + beta + alpha * beta;
//       result_grad_type grad = result_grad_type::Zero();
//
//       grad(0,0) = (1.0/lambda + alpha);
//       grad(0,1) = /* f'(Y)= */ M_PI*alpha * cos(M_PI*p.y());
//
//       grad(1,1) = (1.0/lambda - gamma/(1.0 + gamma));
//
//       grad(2,0) = /* g'(X)= */ M_PI*beta* cos(M_PI*p.x());
//       grad(2,1) = /* h'(Y)= */ 0.0;
//       grad(2,2) = (1.0/lambda + beta);
//
//       return factor*grad;
//    };

   //load and solutions
   auto load = [](const point<T,3>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.0;
      T fz = 0.0;
      return result_type{fx, fy, fz};
   };

   auto solution = [](const point<T,3>& p) -> result_type {
      T fx = -1.0;
      T fy = -1.0;
      T fz = -1.0;

      return result_type{fx,fy,fz};
   };


   auto neumann = [elas_param](const point<T,3>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.0;
      T fz = 0.0;

      return result_type{fx,fy,fz};
   };

   // Define Boundary Conditions
   BoundaryType N1;
   N1.id = 4;//4 cyl
   N1.boundary_type = FREE;

   BoundaryType N2;
   N2.id = 11;
   N2.boundary_type = FREE;

   const std::vector<BoundaryType> boundary_neumann = {N1, N2};

   BoundaryType D1;
   D1.id = 18;
   D1.boundary_type = CLAMPED;

   const std::vector<BoundaryType> boundary_dirichlet = {D1};

   // Solve
   hyperelasticity_solver<Mesh, T, 3, Storage,  point<T, 3> >
      nl(msh, rp, elas_param, boundary_neumann, boundary_dirichlet);

   nl.compute_initial_state();

   if(nl.verbose()){
      std::cout << "Solving the problem ..." << '\n';
   }

   const SolverInfo solve_info = nl.compute(load, solution, neumann);

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

   // PostProcessing
   if(nl.test_convergence()){
      //std::cout << "avrage diameter h: " << average_diameter(msh) << std::endl;
      //std::cout << "l2 error: " << nl.compute_l2_error(solution) << std::endl;

        //std::cout << "Post-processing: " << std::endl;
        //nl.compute_discontinuous_solution("sol_disc3D.msh");
        //nl.compute_deformed("def3D.msh");
   }

   std::cout << "Post-processing: " << std::endl;
   std::string name = "Result_cav_k" + std::to_string(rp.m_cell_degree) + "_l" + std::to_string(rp.m_grad_degree)
   + "_b" + std::to_string(rp.m_beta) + "_s" + std::to_string(rp.m_stab) + "_";
   nl.compute_discontinuous_displacement(name + "sol_disc.msh");
   nl.compute_continuous_displacement(name +"sol_cont.msh");
   nl.compute_deformed(name +"def.msh");
   nl.compute_J_GP(name +"J_GP.msh");
   nl.compute_continuous_J(name +"J_cont.msh");
   nl.compute_discontinuous_J(name +"J_disc.msh");
   nl.compute_discontinuous_Prr(name +"Prr.msh", "Prr");
   nl.compute_discontinuous_Prr(name +"Poo.msh", "Poo");
   nl.compute_discontinuous_VMIS(name +"VM_disc.msh");
   nl.compute_continuous_VMIS(name +"VM_cont.msh");
   nl.compute_VMIS_GP(name +"VM_GP.msh");
}

//Main
int main(int argc, char **argv)
{
    using RealType = double;

    char    *mesh_filename  = nullptr;
    int face_degree, cell_degree, grad_degree, l;

    face_degree = 0; cell_degree = 0; grad_degree = 0;

    ParamRun<RealType> rp;

    ElasticityParameters param = ElasticityParameters();

    param.mu = 0.1;
    param.lambda = 1.0;
    param.type_law = 1;

    //Read Parameters
    int ch;

    while ( (ch = getopt(argc, argv, "g:k:l:p:r:v:w")) != -1 )
    {
        switch(ch)
        {
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

            case 'p':
               rp.m_beta = atof(optarg);
               break;

            case 'r':
               if(!rp.readParameters(optarg))
               exit(1);

               break;

            case 'v':
                rp.m_verbose = true;
                break;

            case 'w': rp.m_stab = false;
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


    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad2$") ))
    {
       std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
       auto msh = disk::load_cartesian_2d_mesh2<RealType>(mesh_filename);
       run_hyperelasticity_solver(msh, rp, param);
       return 0;
    }

    /* Medit 2d*/
    if (std::regex_match(mesh_filename, std::regex(".*\\.medit2d$") ))
    {
       std::cout << "Guessed mesh format: Medit format" << std::endl;
       auto msh = disk::load_medit_2d_mesh<RealType>(mesh_filename);
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

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
    {
       std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
       auto msh = disk::load_fvca6_3d_mesh<RealType>(mesh_filename);
       run_hyperelasticity_solver(msh, rp, param);
       return 0;
   }

}
