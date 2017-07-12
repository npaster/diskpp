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
   size_t  m_degree;        //face degree
   int     m_l;             //cell_degree = face_degree + l 

   size_t  m_n_time_step;   //number of time time_step
   size_t  m_sublevel;      //number od sublevel if there are problems

   bool    m_prediction;    //prediction  u_0(t+dt) = (1 +dt) * u_conv(t)
   bool    m_verbose;       //some printing
   
   size_t  m_iter_max;        //maximun nexton iteration
   T       m_epsilon;         //stop criteria    
   

   
   ParamRun() : m_degree(1), m_l(0), m_n_time_step(1), m_sublevel(1), 
                m_verbose(false), m_prediction(false),
                m_iter_max(10), m_epsilon(T(10E-6)) {}
};
   





