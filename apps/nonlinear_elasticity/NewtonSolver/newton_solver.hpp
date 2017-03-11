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

#include "../../config.h"

#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "hho/hho.hpp"
#include "newton_step.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

struct newton_info
{
    size_t  linear_system_size;
    double  time_assembly, time_solve, time_post;
};


template<typename Mesh>
class NewtonRaphson_solver
{
   typedef Mesh                                        mesh_type;
   typedef typename mesh_type::scalar_type            scalar_type;
   typedef typename mesh_type::cell                   cell_type;
   typedef typename mesh_type::face                   face_type;

   typedef disk::quadrature<mesh_type, cell_type>   cell_quadrature_type;
   typedef disk::quadrature<mesh_type, face_type>   face_quadrature_type;
   typedef disk::scaled_monomial_vector_basis<mesh_type, cell_type>  cell_basis_type;
   typedef disk::scaled_monomial_vector_basis<mesh_type, face_type>  face_basis_type;

   typedef dynamic_matrix<scalar_type>         matrix_dynamic;
   typedef dynamic_vector<scalar_type>         vector_dynamic;

   size_t m_cell_degree, m_face_degree, m_degree;


   const mesh_type& m_msh;

   std::vector<vector_dynamic>                  m_postprocess_data;

   bool m_verbose;
   bool m_convergence;



public:
   NewtonRaphson_solver(const mesh_type& msh, size_t degree, int l = 0)
   : m_msh(msh), m_verbose(false), m_convergence(false)
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


   template<typename LoadIncrement, typename BoundaryConditionFunction>
   newton_info
   compute( const LoadIncrement& lf, const BoundaryConditionFunction& bf,
            const scalar_type epsilon = 1.E-6,
            const std::size_t iter_max = 30, const std::size_t reac_iter=1)
   {
      newton_info ai;
      bzero(&ai, sizeof(ai));
      timecounter tc;

      //initialise the NewtonRaphson_step
      NewtonRaphson_step<Mesh> newton_step(m_msh, m_degree, 0);

      // loop
      std::size_t iter = 0;
      while (iter < iter_max && m_convergence) {
         if(true){
            tc.tic();
            //assemble lhs and rhs
            newton_step.assemble(lf, bf);
            tc.toc();
            ai.time_assembly += tc.to_double();
         }
         // else {
         //    // assemble only rhs
         //    newton_step.assemble_rhs(lf, bf);
         // }
         // solve the global system
         tc.tic();
         newton_step.solve();
         tc.toc();
         ai.time_solve += tc.to_double();
         // update unknowns
         tc.tic();
         //newton_step.update_data();
         tc.toc();
         ai.time_post += tc.to_double();
         // test convergence
         //m_convergence = newton_step.test_convergence(epsilon);
         iter++;
      }

      return ai;
   }

   bool test_convergence() const {return m_convergence;}

   std::vector<scalar_type> result() const
   {
      // export the converged variables
      std::vector<scalar_type> v(4, 0.);
      return v;
   }

};
