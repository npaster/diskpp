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

#pragma once

#include <iostream>

#include <sstream>
#include <string>
#include <math.h>

#include "common/eigen.hpp"
#include "BehaviorLaws/maths_tensor.hpp"
#include <Eigen/LU>

#define _USE_MATH_DEFINES
#include <cmath>

/*
Fichier pour gérer les lois de comportements
-faire une boite qui prends en entrée eps --> behaviorbox --> stress
-penser aux variables internes pour sauvegarder partie elastique
- utilise t on des tenseur pour le module tangent
*/


template<typename scalar_type>
class LawNoLinear
{
   scalar_type m_lambda;


public:
   LawNoLinear()
   :  m_lambda(0.0)
   {}

   LawNoLinear(const scalar_type lambda)
   : m_lambda(lambda)
   {}



   void
   setLambda(const scalar_type lambda)
   {
      m_lambda = lambda;
   }


   scalar_type
   giveLambda() const {return m_lambda;}


   template< int DIM>
   static_matrix<scalar_type, DIM, DIM>
   compute_PK1(const static_matrix<scalar_type, DIM, DIM>& G)
   {
      return  m_lambda * G.trace() * G;
   }


   template<int DIM>
   static_tensor<scalar_type, DIM>
   compute_tangent_moduli_A(const static_matrix<scalar_type, DIM, DIM>& G)
   {
      const static_matrix<scalar_type, DIM, DIM> I2 = static_matrix<scalar_type, DIM, DIM>::Identity();
      const static_tensor<scalar_type, DIM> I4 = compute_IdentityTensor<scalar_type,DIM>();
      const static_tensor<scalar_type, DIM> I2_G = computeKroneckerProduct(I2, G);


      return m_lambda * (I2_G + G.trace() * I4);
   }


   template<int DIM>
   std::pair<static_matrix<scalar_type, DIM, DIM>, static_tensor<scalar_type, DIM> >
   compute_whole_PK1(const static_matrix<scalar_type, DIM, DIM>& G)
   {
      const static_matrix<scalar_type, DIM, DIM> I2 = static_matrix<scalar_type, DIM, DIM>::Identity();
      const static_tensor<scalar_type, DIM> I4 = compute_IdentityTensor<scalar_type,DIM>();
      const static_tensor<scalar_type, DIM> I2_G = computeKroneckerProduct(I2, G);

      const auto PK1 = m_lambda * G.trace() * G;

      const auto A = m_lambda * (I2_G + G.trace() * I4);

      return std::make_pair(PK1, A);
   }

};
