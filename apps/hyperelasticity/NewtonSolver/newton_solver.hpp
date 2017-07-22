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

 // NewtonRaphson_solver

#pragma once

#include <iostream>

#include <sstream>


#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "hho/hho.hpp"
#include "newton_step.hpp"
#include "../ElasticityParameters.hpp"
#include "../BoundaryConditions.hpp"
#include "../Informations.hpp"
#include "../Parameters.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>


template<typename BQData>
class NewtonRaphson_solver_hyperelasticity
{
   typedef typename BQData::mesh_type          mesh_type;
   typedef typename mesh_type::scalar_type     scalar_type;

   typedef dynamic_matrix<scalar_type>         matrix_dynamic;
   typedef dynamic_vector<scalar_type>         vector_dynamic;

   typedef ParamRun<scalar_type>               param_type;
   const BQData&                               m_bqd;

   const mesh_type&                            m_msh;

   std::vector<vector_dynamic>         m_solution_data;
   std::vector<vector_dynamic>         m_solution_cells, m_solution_faces, m_solution_lagr;

   bool m_verbose;
   bool m_convergence;

   BoundaryConditions m_boundary_condition;

   ElasticityParameters m_elas_param;

   const param_type& m_rp;

public:
   NewtonRaphson_solver_hyperelasticity(const mesh_type& msh, const BQData& bqd,
                                        const param_type& rp,
                                        const ElasticityParameters elas_param,
                                        const BoundaryConditions& boundary_conditions)
   : m_msh(msh), m_verbose(false), m_convergence(false), m_elas_param(elas_param), m_bqd(bqd),
     m_rp(rp), m_boundary_condition(boundary_conditions)
   {}

   bool    verbose(void) const     { return m_verbose; }
   void    verbose(bool v)         { m_verbose = v; }


    void
    initialize( const std::vector<vector_dynamic>& initial_solution_cells,
                const std::vector<vector_dynamic>& initial_solution_faces,
                const std::vector<vector_dynamic>& initial_solution_lagr,
                const std::vector<vector_dynamic>& initial_solution)
    {
      m_solution_cells.clear();
      m_solution_cells = initial_solution_cells;
      assert(m_msh.cells_size() == m_solution_cells.size());

      m_solution_faces.clear();
      m_solution_faces = initial_solution_faces;
      assert(m_msh.faces_size() == m_solution_faces.size());

      m_solution_lagr.clear();
      m_solution_lagr = initial_solution_lagr;

      m_solution_data.clear();
      m_solution_data = initial_solution;
      assert(m_msh.cells_size() == m_solution_data.size());

    }

    template<typename LoadIncrement, typename BoundaryConditionFunction, typename NeumannFunction>
    NewtonSolverInfo
   compute( const LoadIncrement& lf, const BoundaryConditionFunction& bf, const NeumannFunction& g)
   {
      NewtonSolverInfo ni;
      timecounter tc;
      tc.tic();

      const size_t iter_max = m_rp.m_iter_max;
      const scalar_type epsilon = m_rp.m_epsilon;
      bool auricchio = false;
      scalar_type error;

      //initialise the NewtonRaphson_step
      NewtonRaphson_step_hyperelasticity<BQData>
         newton_step(m_msh, m_bqd, m_rp, m_elas_param, m_boundary_condition);

      newton_step.initialize(m_solution_cells, m_solution_faces, m_solution_lagr, m_solution_data);
      newton_step.verbose(m_verbose);

      m_convergence = false;

      size_t nb_negative_ev_init = 0;
      // loop
      std::size_t iter = 0;
      while (iter < iter_max && !m_convergence) {

          //assemble lhs and rhs
          AssemblyInfo assembly_info;
          try {
             assembly_info = newton_step.assemble(lf, bf, g);
          }
          catch(const std::invalid_argument& ia){
                std::cerr << "Invalid argument: " << ia.what() << std::endl;
                m_convergence = false;
                tc.toc();
                ni.m_time_newton = tc.to_double();
                return ni;
          }

          ni.updateAssemblyInfo( assembly_info);
         // test convergence
         m_convergence = newton_step.test_convergence(epsilon, iter, error);

         if(auricchio){
            size_t nb_negative_ev = newton_step.test_aurrichio();
            if(iter == 0)
               nb_negative_ev_init = nb_negative_ev;
            else if(nb_negative_ev > nb_negative_ev_init)
               std::cout << "Test Aurricchio: we loos the coercivite of D2L " << nb_negative_ev << " > " << nb_negative_ev_init << std::endl;
         }

         if(iter < (iter_max-1) && !m_convergence){
            // solve the global system
            SolveInfo solve_info = newton_step.solve();
            ni.updateSolveInfo(solve_info);
            // update unknowns
            PostprocessInfo post_info = newton_step.postprocess(lf);
            ni.updatePostProcessInfo(post_info);
            newton_step.update_solution();

            if(m_rp.m_compute_energy and m_verbose){
               std::array<scalar_type,2> energy = newton_step.compute_energy();
               std::cout << "Compute elastic energy:" << std::endl;
               std::cout << " - Elastic energy: " << energy[0]  << std::endl;
               std::cout << " - Stabilisation energy: " << energy[1]  << std::endl;
            }
         }
         iter++;
      }

//       if(!m_convergence)
//          m_convergence = newton_step.test_convergence(1.E-6, iter_max);

      if(m_convergence)
         newton_step.save_solutions(m_solution_cells, m_solution_faces, m_solution_lagr, m_solution_data);

      tc.toc();
      ni.m_time_newton = tc.to_double();
      return ni;
   }

   bool test_convergence() const {return m_convergence;}

   void
    save_solutions( std::vector<vector_dynamic>& solution_cells,
                    std::vector<vector_dynamic>& solution_faces,
                    std::vector<vector_dynamic>& solution_lagr,
                    std::vector<vector_dynamic>& solution)
    {
      solution_cells.clear();
      solution_cells = m_solution_cells;
      assert(m_solution_cells.size() == solution_cells.size());

      solution_faces.clear();
      solution_faces = m_solution_faces;
      assert(m_solution_faces.size() == solution_faces.size());

      solution_lagr.clear();
      solution_lagr = m_solution_lagr;
      assert(m_solution_lagr.size() == solution_lagr.size());

      solution.clear();
      solution = m_solution_data;
      assert(m_solution_data.size() == solution.size());
    }

};
