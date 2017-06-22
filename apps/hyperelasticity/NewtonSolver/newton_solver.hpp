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

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

struct newton_info
{
    size_t  linear_system_size;
    double  time_assembly, time_solve, time_post;
    double  time_gradrec, time_statcond, time_stab, time_elem, time_law;
};


template<typename Mesh>
class NewtonRaphson_solver_hyperelasticity
{
   typedef Mesh                                        mesh_type;
   typedef typename mesh_type::scalar_type            scalar_type;



   typedef dynamic_matrix<scalar_type>         matrix_dynamic;
   typedef dynamic_vector<scalar_type>         vector_dynamic;

   size_t m_cell_degree, m_face_degree, m_degree;


   const mesh_type& m_msh;

   std::vector<vector_dynamic>         m_solution_data;
   std::vector<vector_dynamic>         m_solution_cells, m_solution_faces, m_solution_lagr;

   bool m_verbose;
   bool m_convergence;

   ElasticityParameters m_elas_param;



public:
   NewtonRaphson_solver_hyperelasticity(const mesh_type& msh, const size_t degree, const ElasticityParameters elas_param, int l = 0)
   : m_msh(msh), m_verbose(false), m_convergence(false), m_elas_param(elas_param)
   {
      if ( l < -1 or l > 1)
      {
         std::cout << "'l' should be -1, 0 or 1. Reverting to 0." << std::endl;
         l = 0;
      }

      if (degree == 0 && l == -1)
      {
         std::cout << "'l' should be 0 or 1. Reverting to 0." << std::endl;
         l = 0;
      }

      m_cell_degree = degree + l;
      m_face_degree = degree;
      m_degree = degree;

   }

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
      assert((m_msh.boundary_faces_size() - number_of_neumann_faces(m_msh, boundary_neumann))
            == m_solution_lagr.size());

      m_solution_data.clear();
      m_solution_data = initial_solution;
      assert(m_msh.cells_size() == m_solution_data.size());

    }

    template<typename LoadIncrement, typename BoundaryConditionFunction, typename NeumannFunction>
   newton_info
   compute( const LoadIncrement& lf, const BoundaryConditionFunction& bf, const NeumannFunction& g,
            const std::vector<size_t>& boundary_neumann,
            const scalar_type epsilon = 1.E-7,
            const std::size_t iter_max = 20)
   {
      newton_info ai;
      bzero(&ai, sizeof(ai));
      timecounter tc;

      //initialise the NewtonRaphson_step
      NewtonRaphson_step_hyperelasticity<Mesh> newton_step(m_msh, m_degree, m_elas_param, 0);

      newton_step.initialize(m_solution_cells, m_solution_faces,
                             m_solution_lagr, m_solution_data);
      newton_step.verbose(m_verbose);

      m_convergence = false;
      // loop
      std::size_t iter = 0;
      while (iter < iter_max && !m_convergence) {
          tc.tic();
          //assemble lhs and rhs
          auto time_assembly = newton_step.assemble(lf, bf, g, boundary_neumann);
          tc.toc();
          ai.time_assembly += tc.to_double();

          ai.time_gradrec +=  time_assembly.time_gradrec;
          ai.time_stab +=  time_assembly.time_stab;
          ai.time_statcond +=  time_assembly.time_statcond;
          ai.time_elem +=  time_assembly.time_elem;
          ai.time_law +=  time_assembly.time_law;
         // test convergence
         m_convergence = newton_step.test_convergence(epsilon, iter);

         if(iter < (iter_max-1) && !m_convergence){
            // solve the global system
            tc.tic();
            solver_info solve_info = newton_step.solve();
            tc.toc();
            ai.time_solve += solve_info.time_solver;
            // update unknowns
            tc.tic();
            newton_step.postprocess(lf, g, boundary_neumann);
            newton_step.update_solution();
            tc.toc();
            ai.time_post += tc.to_double();
         }
         iter++;
      }

      if(!m_convergence)
         m_convergence = newton_step.test_convergence(1.E-6, iter_max);

        newton_step.save_solutions(m_solution_cells, m_solution_faces,
                             m_solution_lagr, m_solution_data);

      return ai;
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
