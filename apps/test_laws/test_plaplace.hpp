/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
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

#include "common/eigen.hpp"
#include "BehaviorLaws/behaviorlaws.hpp"

#pragma once


template< typename T>
void
test_plaplace_2d(const T prec)
{
   for( size_t p = 2; p < 5; p++)
   {
      static_vector<T, 2> G;
      G(0) = T{1.0}; G(1) = std::sqrt(3.0);

      T G_norm2 = std::pow(2.0, p-2.0);
      T G_norm4 = std::pow(2.0, p-4.0);

      static_matrix<T, 2, 2> A;
      A(0,0) = T{1.0}; A(1,0) = std::sqrt(3.0);
      A(0,1) = std::sqrt(3.0); A(1,1) = T{3.0};

      A *= (p-2) * G_norm4;

      A(0,0) += G_norm2; A(1,1) += G_norm2;

      pLaplaceLaw<T> law(p);

      auto Ac = law.compute_tangent_moduli(G);

      T error = (A - Ac).norm();

      if(error > prec)
         std::cout << "pLaplace 2d, error A: " << error << std::endl;
   }

}


template< typename T>
void
test_plaplace_3d(const T prec)
{
   for( size_t p = 2; p < 5; p++)
   {
      static_vector<T, 3> G;
      G(0) = T{1.0}; G(1) = std::sqrt(3.0); G(2) = std::sqrt(5.0);

      T G_norm2 = std::pow(3.0, p-2.0);
      T G_norm4 = std::pow(3.0, p-4.0);

      static_matrix<T, 3, 3> A;
      A(0,0) = T{1.0}; A(0,1) = std::sqrt(3.0); A(0,2) = std::sqrt(5.0);
      A(1,0) = std::sqrt(3.0); A(1,1) = T{3.0}; A(1,2) = std::sqrt(15.0);
      A(2,0) = std::sqrt(5.0); A(2,1) = std::sqrt(15.0); A(2,2) = 5.0;

      A *= (p-2) * G_norm4;

      A(0,0) += G_norm2; A(1,1) += G_norm2; A(2,2) += G_norm2;

      pLaplaceLaw<T> law(p);

      auto Ac = law.compute_tangent_moduli(G);

      T error = (A - Ac).norm();

      if(error > prec)
         std::cout << "pLaplace 3d, error A: " << error << std::endl;
   }

}
