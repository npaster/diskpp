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
 * cite
*/

#pragma once


template< typename T>
class ParamRun
{
public:
   size_t  m_face_degree;     //face degree
   size_t  m_cell_degree;    //cell_degree = face_degree + l
   size_t  m_grad_degree;     // grad degree
   int     m_l;

   size_t  m_n_time_step;   //number of time time_step
   size_t  m_sublevel;      //number od sublevel if there are problems

   bool    m_stab;          //stabilization yes or no
   bool    m_adapt_coeff;   //adapts automatically the stabilisation coefficient
   bool    m_adapt_stab;    //use the adpatative stabilization
   bool    m_verbose;       //some printing
   bool    m_compute_energy; // to compute intere energy

   bool    m_init;          // a first time step
   T       m_tinit;         // time of the first time step

   size_t  m_iter_max;        //maximun nexton iteration
   T       m_epsilon;         //stop criteria

   T       m_beta_min;      // minium of stabilization constant
   T       m_beta_max;      // maximum of stabilization constant


   ParamRun() : m_face_degree(1), m_cell_degree(1), m_grad_degree(1), m_l(0),
                m_n_time_step(1), m_sublevel(1), m_stab(true),
                m_verbose(false), m_adapt_coeff(false), m_adapt_stab(false),
                m_compute_energy(false), m_init(false), m_tinit(0),
                m_iter_max(10), m_epsilon(T(10E-6)),
                m_beta_min(T(0)), m_beta_max(T(0)) {}
};
