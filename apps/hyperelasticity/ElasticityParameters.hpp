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



 struct ElasticityParameters
   {
      ElasticityParameters()
      {
         double lambda = 1.0;
         double mu = 1.0;
         double tau = 1.0;
         bool adaptative_stab = false;
         size_t type_law = 1;
      }
         double lambda;
         double mu;
         double tau; // stabilisation parameter
         bool adaptative_stab;
         size_t type_law;
   };
   
   
 enum BoundaryType : size_t
 {
    INTERNAL = 0,
    DIRICHLET = 1,
    NEUMANN = 2,
 };




