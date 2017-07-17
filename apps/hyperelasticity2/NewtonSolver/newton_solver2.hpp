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
#include "newton_step2.hpp"
#include "../ElasticityParameters.hpp"
#include "../BoundaryConditions.hpp"
#include "../Informations.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>


template<typename BQData>
class NewtonRaphson_solver_hyperelasticity2
{
   typedef typename BQData::mesh_type          mesh_type;
   typedef typename mesh_type::scalar_type     scalar_type;
   
   typedef dynamic_matrix<scalar_type>         matrix_dynamic;
   typedef dynamic_vector<scalar_type>         vector_dynamic;
   
   const BQData&                               m_bqd;

   const mesh_type&                            m_msh;

   std::vector<vector_dynamic>         m_solution_data;
   std::vector<vector_dynamic>         m_solution_cells, m_solution_faces, m_solution_lagr;

   bool m_verbose;
   bool m_convergence;

   ElasticityParameters m_elas_param;


public:
   NewtonRaphson_solver_hyperelasticity2(const mesh_type& msh, const BQData& bqd, const ElasticityParameters elas_param)
   : m_msh(msh), m_verbose(false), m_convergence(false), m_elas_param(elas_param), m_bqd(bqd)
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
    
    void
    prediction(const scalar_type delta_t)
    {
       scalar_type incr = scalar_type(1.0 + delta_t);
       
       for(size_t i = 0; i < m_solution_cells.size(); i++)
          m_solution_cells[i] *= incr;
       
       for(size_t i = 0; i < m_solution_faces.size(); i++)
          m_solution_faces[i] *= incr;
       
       for(size_t i = 0; i < m_solution_lagr.size(); i++)
          m_solution_lagr[i] *= incr;
       
       for(size_t i = 0; i < m_solution_data.size(); i++)
          m_solution_data[i] *= incr;
       
    }

    template<typename LoadIncrement, typename BoundaryConditionFunction, typename NeumannFunction>
    NewtonSolverInfo
   compute( const LoadIncrement& lf, const BoundaryConditionFunction& bf, const NeumannFunction& g,
            const std::vector<size_t>& boundary_neumann, const std::vector<BoundaryConditions>& boundary_dirichlet,
            const scalar_type epsilon = 1.E-6,
            const std::size_t iter_max = 10)
   {
      NewtonSolverInfo ni;
      timecounter tc;
      tc.tic();
      
      bool auricchio = false;

      //initialise the NewtonRaphson_step
      NewtonRaphson_step_hyperelasticity2<BQData> newton_step(m_msh, m_bqd, m_elas_param);

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
             assembly_info = newton_step.assemble(lf, bf, g, boundary_neumann, boundary_dirichlet);
          }
          catch(const std::invalid_argument& ia){
                std::cerr << "Invalid argument: " << ia.what() << std::endl;
                m_convergence = false;
                tc.toc();
                ni.m_time_newton = tc.to_double();
                return ni;
          }
          
          //je pense que l'on peut le supprimer
          if(auricchio && iter == 0){
             nb_negative_ev_init = newton_step.test_aurrichio();
          }

          ni.updateAssemblyInfo( assembly_info);
         // test convergence
          scalar_type error(0.0);
         m_convergence = newton_step.test_convergence(epsilon, iter, error);


         if(iter < (iter_max-1) && !m_convergence){
            // solve the global system
            SolveInfo solve_info = newton_step.solve();
            ni.updateSolveInfo(solve_info);
            // update unknowns
            PostprocessInfo post_info = newton_step.postprocess(lf, g, boundary_neumann, boundary_dirichlet);
            ni.updatePostProcessInfo(post_info);

            newton_step.update_solution();
         }
         iter++;
      }

//       if(!m_convergence)
//          m_convergence = newton_step.test_convergence(1.E-6, iter_max);

      if(m_convergence)
         newton_step.save_solutions(m_solution_cells, m_solution_faces, m_solution_lagr, m_solution_data);
      
      
      if(auricchio && m_convergence){
         size_t nb_negative_ev = newton_step.test_aurrichio();
         if(nb_negative_ev > nb_negative_ev_init)
            std::cout << "Test Aurricchio: we loos the coercivite of D2L " << nb_negative_ev << " > " << nb_negative_ev_init << std::endl;
      }

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