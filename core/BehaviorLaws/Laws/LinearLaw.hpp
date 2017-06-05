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

#include "common/eigen.hpp"
#include "BehaviorLaws/maths_tensor.hpp"


#define _USE_MATH_DEFINES
#include <cmath>

/*
Fichier pour gérer les lois de comportements
-faire une boite qui prends en entrée eps --> behaviorbox --> stress
-penser aux variables internes pour sauvegarder partie elastique
- utilise t on des tenseur pour le module tangent
*/

/* Material: LinearLaw
 * Stress :  S(G) = lambda * G
 * Module :  A(G) = lambda * I_4
 */


template<typename scalar_type>
class LinearLaw
{
   typedef scalar_type   tensor_1d;
   typedef static_tensor<scalar_type, 2>   tensor_2d;
   typedef static_tensor<scalar_type, 3>   tensor_3d;
   scalar_type m_lambda;

   tensor_1d  linear_tensor_1d;
   tensor_2d  linear_tensor_2d;
   tensor_3d  linear_tensor_3d;

   void
   init_linear_tensor()
   {
      linear_tensor_1d = m_lambda;
      linear_tensor_2d = m_lambda * compute_IdentityTensor<scalar_type,2>();
      linear_tensor_3d = m_lambda * compute_IdentityTensor<scalar_type,3>();

   }

public:
   LinearLaw()
   : m_lambda(0.0)
   {
      linear_tensor_1d = 0.0;
      linear_tensor_2d.setZero();
      linear_tensor_3d.setZero();
   }

   LinearLaw(const scalar_type lambda)
   : m_lambda(lambda)
   {
      init_linear_tensor();
   }

   void
   setLambda(const scalar_type lambda)
   {
      m_lambda = lambda;
      init_linear_tensor();
   }

   scalar_type
   giveLambda() const {return m_lambda;}


   template<int DIM>
   static_matrix<scalar_type, DIM, DIM>
   compute_P(const static_matrix<scalar_type, DIM, DIM>& Gradient)
   {
      return m_lambda * Gradient;
   }
   
   //1D case
   scalar_type
   compute_P(const scalar_type& Gradient)
   {
      return m_lambda * Gradient;
   }

   tensor_1d
   compute_tangent_moduli(const static_matrix<scalar_type,1,1>& Gradient)
   {
      return linear_tensor_1d;
   }
   
   tensor_2d
   compute_tangent_moduli(const static_matrix<scalar_type,2,2>& Gradient)
   {
      return linear_tensor_2d;
   }

   tensor_3d
   compute_tangent_moduli(const static_matrix<scalar_type,3,3>& Gradient)
   {
      return linear_tensor_3d;
   }

   std::pair<scalar_type, scalar_type>
   compute_whole(const scalar_type& Gradient)
   {
      scalar_type P = compute_P(Gradient);
      
      return std::make_pair(P, linear_tensor_1d);
   }

   std::pair<static_matrix<scalar_type, 2, 2>, static_tensor<scalar_type, 2> >
   compute_whole(const static_matrix<scalar_type, 2, 2>& Gradient)
   {
      static_matrix<scalar_type, 2, 2> P = compute_P(Gradient);

      return std::make_pair(P, linear_tensor_2d);
   }

   std::pair<static_matrix<scalar_type, 3, 3>, static_tensor<scalar_type, 3> >
   compute_whole(const static_matrix<scalar_type, 3, 3>& Gradient)
   {
      static_matrix<scalar_type, 3, 3> P = compute_P(Gradient);

      return std::make_pair(P, linear_tensor_3d);
   }

};
