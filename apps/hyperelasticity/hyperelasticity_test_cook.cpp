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

struct resultat_type
{
   size_t  degree;
   size_t nb_dof;
   double h;
   double ux;
   double uy;
   std::string name;
};

template<template<typename, size_t , typename> class Mesh,
         typename T, typename Storage>
resultat_type
run_cook_solver(const Mesh<T, 2, Storage>& msh, const ParamRun<T>& rp, const ElasticityParameters& elas_param,
                 const std::string& name)
{
   typedef Mesh<T, 2, Storage> mesh_type;
   typedef static_vector<T, 2> result_type;
   typedef static_matrix<T, 2, 2> result_grad_type;

   auto load = [](const point<T,2>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.0;;
      return result_type{fx,fy};
   };

   auto solution = [](const point<T,2>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.0;

      return result_type{fx,fy};
   };


   auto neumann = [](const point<T,2>& p) -> result_type {
      T fx = 0.0;
      T fy = 0.1;

      return result_type{fx,fy};
   };

   //NEUMANN CONDITION

   BoundaryType N1;
   N1.id = 3;
   N1.boundary_type = FREE;

   BoundaryType N2;
   N2.id = 8;
   N2.boundary_type = FREE;

   BoundaryType N3;
   N3.id = 6;
   N3.boundary_type = NEUMANN;

   std::vector<BoundaryType> boundary_neumann = {N1, N2, N3};

   std::vector<BoundaryType> boundary_dirichlet = {};

   hyperelasticity_solver<Mesh, T, 2, Storage,  point<T, 2> >
      nl(msh, rp, elas_param, boundary_neumann, boundary_dirichlet);


   nl.compute_initial_state();

   if(nl.verbose()){
      std::cout << "Solving the problem ..."  << '\n';
   }

   SolverInfo solve_info = nl.compute(load, solution, neumann);

   if(nl.verbose()){
      std::cout << " " << std::endl;
      std::cout << "------------------------------------------------------- " << std::endl;
      std::cout << "Summaring: " << std::endl;
      std::cout << "Total time to solve the problem: " << solve_info.m_time_solver << " sec" << std::endl;
      std::cout << "**** Assembly time: " << solve_info.m_newton_info.m_assembly_info.m_time_assembly << " sec" << std::endl;
      std::cout << "****** Gradient reconstruction: " << solve_info.m_newton_info.m_assembly_info.m_time_gradrec << " sec" << std::endl;
      std::cout << "****** Elementary computation: " << solve_info.m_newton_info.m_assembly_info.m_time_elem << " sec" << std::endl;
      std::cout << "       *** Behavior computation: " << solve_info.m_newton_info.m_assembly_info.m_time_law << " sec" << std::endl;
      std::cout << "****** Static condensation: " << solve_info.m_newton_info.m_assembly_info.m_time_statcond << " sec" << std::endl;
      std::cout << "**** Solver time: " << solve_info.m_newton_info.m_solve_info.m_time_solve << " sec" << std::endl;
      std::cout << "**** Postprocess time: " << solve_info.m_newton_info.m_postprocess_info.m_time_post << " sec" << std::endl;
      std::cout << "------------------------------------------------------- " << std::endl;
      std::cout << " " << std::endl;
   }

   resultat_type resultat;
   resultat.name = name;
   resultat.h = average_diameter(msh);
   resultat.degree = rp.m_degree;
   resultat.nb_dof = nl.getDofs();
   resultat.ux = 0.0;
   resultat.uy = 0.0;

   if(nl.test_convergence()){
        std::cout << "Post-processing: " << std::endl;
        nl.compute_discontinuous_solution(name + "_k_" + std::to_string(rp.m_degree) + "_sol_disc2D.msh");
        nl.compute_conforme_solution(name + "_k_" + std::to_string(rp.m_degree) + "_sol_conf2D.msh");
        nl.plot_J_at_gausspoint(name + "_k_" + std::to_string(rp.m_degree) + "_J_gp2D.msh");

        auto depl = nl.displacement_node(2);

        resultat.ux = depl[0];
        resultat.uy = depl[1];
   }

   return resultat;
}

void
printResults(const std::vector<resultat_type>& resultat)
{
   if(resultat.size() > 0){
      std::ios::fmtflags f( std::cout.flags() );
      std::cout.precision(3);
      std::cout.setf(std::iostream::scientific, std::iostream::floatfield);

      std::cout << "Cook's membrane test for k = " << resultat[0].degree << std::endl;
      std::cout << "-----------------------------------------------------------------------------" << std::endl;
      std::cout << "|  Name mesh   | Size mesh   | Displacement A | Displacement A |    Total   |" << std::endl;
      std::cout << "|              |    h        |       ux       |       uy       | faces DOF  |" << std::endl;
      std::cout << "-----------------------------------------------------------------------------" << std::endl;

      for(size_t i = 0; i < resultat.size(); i++){
         std::string s_dof = " " + std::to_string(resultat[i].nb_dof) + "                  ";
         s_dof.resize(10);

         std::string s_name = " " + resultat[i].name + "                  ";
         s_name.resize(14);

         std::cout << "|" <<  s_name << "|  " << resultat[i].h << "  |   " << resultat[i].ux << "   |    " <<
         resultat[i].uy <<  "   | "  << s_dof <<  " |" << std::endl;
      }

      std::cout << "-----------------------------------------------------------------------------" << std::endl;
      std::cout << "  " <<std::endl;
      std::cout.flags( f );
   }
   else
      std::cout << "The file error is empty" << std::endl;
}

template< typename T>
void test_cook(const ParamRun<T>& rp, const ElasticityParameters& elas_param)
{
   size_t runs = 5;

   std::vector<std::string> paths;
   paths.push_back("../diskpp/meshes/Tests/Cook/Cook-2x2.medit2d");
   paths.push_back("../diskpp/meshes/Tests/Cook/Cook-4x4.medit2d");
   paths.push_back("../diskpp/meshes/Tests/Cook/Cook-8x8.medit2d");
   paths.push_back("../diskpp/meshes/Tests/Cook/Cook-16x16.medit2d");
   paths.push_back("../diskpp/meshes/Tests/Cook/Cook-32x32.medit2d");
   paths.push_back("../diskpp/meshes/Tests/Cook/Cook-64x64.medit2d");
   paths.push_back("../diskpp/meshes/Tests/Cook/Cook-128x128.medit2d");
   paths.push_back("../diskpp/meshes/Tests/Cook/Cook-256x256.medit2d");

   std::vector<std::string> names;
   names.push_back("Cook-2x2");
   names.push_back("Cook-4x4");
   names.push_back("Cook-8x8");
   names.push_back("Cook-16x16");
   names.push_back("Cook-32x32");
   names.push_back("Cook-64x64");
   names.push_back("Cook-128x128");
   names.push_back("Cook-256x256");



   std::vector<resultat_type> results;

   for(int i = 0; i < runs; i++){
      auto msh = disk::load_medit_2d_mesh<T>(paths[i].c_str());
      auto resultat = run_cook_solver(msh, rp, elas_param, names[i]);
      results.push_back(resultat);
   }
   printResults(results);
}


int main(int argc, char **argv)
{
    using RealType = double;

    char    *mesh_filename  = nullptr;
    int degree, n_time_step, l;

    ParamRun<RealType> rp;
    rp.m_sublevel = 8;
    rp.m_verbose = true;

    ElasticityParameters param = ElasticityParameters();

    param.mu = 0.4225;
    param.lambda = 915.4;
    param.type_law = 6;
    param.tau = 100.0;
    param.adaptative_stab = false;

    int ch;

    while ( (ch = getopt(argc, argv, "e:i:k:l:n:p:v")) != -1 )
    {
        switch(ch)
        {
            case 'e': rp.m_epsilon = atof(optarg); break;
            case 'i': rp.m_iter_max = atoi(optarg); break;
            case 'k':
                degree = atoi(optarg);
                if (degree < 0)
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

            case 'n':
               n_time_step = atoi(optarg);
               if (n_time_step == 0)
                {
                    std::cout << "Number of time step must be positive. Falling back to 1." << std::endl;
                    n_time_step = 1;
                }
                rp.m_n_time_step = n_time_step;
                break;


            case 'v':
                rp.m_verbose = true;
                break;


            case 'p':
               param.tau = atof(optarg);
               break;


            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;


    test_cook<RealType>(rp, param);
    return 0;

}
