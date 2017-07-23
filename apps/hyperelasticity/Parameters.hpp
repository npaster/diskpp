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

#include <fstream>
#include <iostream>
#include <string>

enum StabilizationType : size_t
{
   L2 = 0,
   PIKF = 1,
   HHO = 2,
   NOTHING = 3
};

template< typename T>
class ParamRun
{
public:
   size_t  m_face_degree;     //face degree
   size_t  m_cell_degree;    //cell_degree = face_degree + l
   size_t  m_grad_degree;     // grad degree
   int     m_l;

   std::vector<std::pair<T, size_t>>  m_time_step;   //number of time time_step
   size_t  m_sublevel;      //number od sublevel if there are problems

   bool    m_stab;          //stabilization yes or no
   bool    m_adapt_coeff;   //adapts automatically the stabilisation coefficient
   bool    m_adapt_stab;    //use the adpatative stabilization
   bool    m_verbose;       //some printing
   bool    m_compute_energy; // to compute intere energy

   size_t  m_stab_type;      //type of stabilization
   size_t  m_iter_max;        //maximun nexton iteration
   T       m_epsilon;         //stop criteria

   T       m_beta_min;      // minium of stabilization constant
   T       m_beta_max;      // maximum of stabilization constant


   ParamRun() : m_face_degree(1), m_cell_degree(1), m_grad_degree(1), m_l(0),
                m_sublevel(1), m_stab(true), m_stab_type(HHO),
                m_verbose(false), m_adapt_coeff(false), m_adapt_stab(false),
                m_compute_energy(false),
                m_iter_max(10), m_epsilon(T(1E-6)),
                m_beta_min(T(0)), m_beta_max(T(0))
                {
                   m_time_step.push_back(std::make_pair(1.0, 10));
                }


   bool readParameters(const std::string& filename)
   {
      std::ifstream   ifs(filename);
      std::string     keyword;

      if (!ifs.is_open())
      {
         std::cout << "Error opening " << filename << std::endl;
         return false;
      }

      ifs >> keyword;
      if ( keyword != "BeginParameters" )
      {
         std::cout << "Expected keyword \"BeginParameters\"" << std::endl;
         return false;
      }


      ifs >> keyword;
      while( keyword != "EndParameters")
      {
         if ( keyword == "FaceDegree" )
         {
            ifs >> m_face_degree;
         }
         else if ( keyword == "CellDegree" )
         {
            ifs >> m_cell_degree;
         }
         else if ( keyword == "GradDegree" )
         {
            ifs >> m_grad_degree;
         }
         else if ( keyword == "Sublevel" )
         {
            ifs >> m_sublevel;
         }
         else if ( keyword == "TimeStep" )
         {
            size_t n_time_step(0);
            ifs >> n_time_step;

            m_time_step.clear();
            m_time_step.reserve(n_time_step);
            for (size_t i = 0; i < n_time_step; i++) {
               T  time(0.0);
               size_t   time_step(0);
               ifs >> time >> time_step;
               m_time_step.push_back(std::make_pair(time, time_step));
            }
         }
         else if ( keyword == "Stabilisation" )
         {
            std::string logical;
            ifs >> logical;
            if(logical == "true")
               m_stab = true;
            else{
               m_stab = false;
               m_stab_type = NOTHING;
            }
         }
         else if ( keyword == "TypeStabiliszation" )
         {
            std::string type;
            ifs >> type;

            if(type == "L2")
               m_stab_type = L2;
            else if(type == "PIKF")
               m_stab_type = PIKF;
            else if(type == "HHO")
               m_stab_type = HHO;
            else if(type == "NOTHING")
               m_stab_type = NOTHING;
         }
         else if ( keyword == "AdaptativeCoefficient" )
         {
            std::string logical;
            ifs >> logical;

            if(logical == "true")
               m_adapt_coeff = true;
            else
               m_adapt_coeff = false;
         }
         else if ( keyword == "AdaptativeStabilization" )
         {
            std::string logical;
            ifs >> logical;

            if(logical == "true")
               m_adapt_stab = true;
            else
               m_adapt_stab = false;
         }
         else if ( keyword == "BetaMax" )
         {
            ifs >> m_beta_max;
         }
         else if ( keyword == "BetaMin" )
         {
            ifs >> m_beta_min;
         }
         else if ( keyword == "Verbose" )
         {
            std::string logical;
            ifs >> logical;
            if(logical == "true")
               m_verbose = true;
            else
               m_verbose = false;
         }
         else if ( keyword == "ComputeEnergy" )
         {
            std::string logical;
            ifs >> logical;

            if(logical == "true")
               m_compute_energy = true;
            else
               m_compute_energy = false;
         }
         else if ( keyword == "IterMax" )
         {
            ifs >> m_iter_max;
         }
         else if ( keyword == "Epsilon" )
         {
            ifs >> m_epsilon;
         }
         else
         {
            std::cout << "Error parsing Parameters file" << std::endl;
            return false;
         }

         ifs >> keyword;
      }

      ifs.close();
      return true;
   }
};
