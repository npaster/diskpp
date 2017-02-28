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

#include <iostream>

#include <sstream>

#include "../../config.h"

#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "hho/hho.hpp"
##include "newton_step.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

struct assembly_info
{
    size_t  linear_system_size;
    double  time_gradrec, time_statcond, time_stab;
};

struct solver_info
{
    double  time_solver;
};

struct postprocess_info
{
    double  time_postprocess;
};


template<typename Mesh>
class NewtonRaphson_solver
{
    typedef Mesh                                       mesh_type;
    typedef typename mesh_type::scalar_type            scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    size_t m_cell_degree, m_face_degree;

    bqdata_type     m_bqd;

    const mesh_type& m_msh;

    std::vector<dynamic_vector<scalar_type>>        m_postprocess_data;

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

        m_bqd = bqdata_type(m_cell_degree, m_face_degree);
    }

    bool    verbose(void) const     { return m_verbose; }
    void    verbose(bool v)         { m_verbose = v; }


    template<typename LoadIncrement, typename BoundaryConditionFunction>
    void
    compute(const LoadIncrement& lf, const BoundaryConditionFunction& bf,
            const scalar_type epsilon = 1.E-6,
            const std::size_t iter_max = 30, const std::size_t reac_iter=1)
   {
      //initialise the NewtonRaphson_step
      NewtonRaphson_step newton_step(m_msh, m_face_degree, 0);

      // loop
      std::size_t iter = 0;
      while (iter < iter_max && m_convergence) {
         if(true){
            //assemble lhs and rhs
            newton_step.assemble(lf, bf);
         }
         // else {
         //    // assemble only rhs
         //    newton_step.assemble_rhs(lf, bf);
         // }
         // solve the global system
         newton_step.solve();
         // update unknowns
         newton_step.update_data();
         // test convergence
         m_convergence = newton_step.test_convergence(epsilon);
         iter++;
      }
   }

   bool test_convergence() const {return m_convergence;}

   std::vector<scalar_type> result() const
   {
      // export the converged variables
      std::vector<scalar_type> v(4, 0.);
      return v;
   }

};
