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
#include "BehaviorLaws/maths_jacobian_derivate.hpp"


#define _USE_MATH_DEFINES
#include <cmath>

/*
Fichier pour gérer les lois de comportements
-faire une boite qui prends en entrée eps --> behaviorbox --> stress
-penser aux variables internes pour sauvegarder partie elastique
- utilise t on des tenseur pour le module tangent
*/

/* Material: LaplaceLaw
 * Stress :  S(G) = lambda(G) * G
 * Module :  A(G) = d lambda / dG \times G  + lambda(G) * I_4
 */

/*
 * Case 1 : lambda(G) = 1; d lambda / dG = 0;
 * Case 2 : lambda(G) = det(G); d lambda / dG = det(G) * G^{-T};
 */


template<typename scalar_type>
class LaplaceLaw
{
   size_t m_type;
   
   const size_t maxtype = 3;
   
   template< int DIM>
   scalar_type
   computeLambda(const static_matrix<scalar_type, DIM, DIM>& G)
   {
      if(m_type == 1)
         return 1.0;
      else if(m_type == 2)
         return G.determinant();
      else if(m_type == 3)
         return G.trace();
      else
         assert(false);
   }
   
   template< int DIM>
   static_matrix<scalar_type, DIM, DIM>
   computeLambdaprime(const static_matrix<scalar_type, DIM, DIM>& G)
   {
      if(m_type == 1)
         return static_matrix<scalar_type, DIM, DIM>::Zero();
      else if(m_type == 2)
         return computeJacobianFirstDerivate(G);
      else if(m_type == 3)
         return static_matrix<scalar_type, DIM, DIM>::Identity();
      else
         assert(false);
   }

public:
   LaplaceLaw()
   : m_type(1)
   {
      if(m_type <= 0 || m_type > maxtype)
      {
         std::cout << "Unknown option for Laplacian Law" << '\n';
         std::cout << "We use lambda(G) = 1" << '\n';
         
         m_type = 1;
      }
   }
   
   LaplaceLaw(const size_t type)
   : m_type(type)
   {
      if(m_type <= 0 || m_type > maxtype)
      {
         std::cout << "Unknown option for Laplacian Law" << '\n';
         std::cout << "We use lambda(G) = 1" << '\n';
         m_type = 1;    
      }
   }
   
   void
   setType(const size_t type)
   {
      m_type = type;
   }

   template< int DIM>
   static_matrix<scalar_type, DIM, DIM>
   compute_P(const static_matrix<scalar_type, DIM, DIM>& G)
   {
      return  computeLambda(G) * G;
   }
   
   
   
   template<int DIM>
   static_tensor<scalar_type, DIM>
   compute_tangent_moduli(const static_matrix<scalar_type, DIM, DIM>& G)
   {
      auto A1 = computeKroneckerProduct(G, computeLambdaprime(G));
      auto A2 = computeLambda(G) * compute_IdentityTensor<scalar_type,DIM>();

      return A1 + A2;
   }
   
   
   template<int DIM>
   std::pair<static_matrix<scalar_type, DIM, DIM>, static_tensor<scalar_type, DIM> >
   compute_whole(const static_matrix<scalar_type, DIM, DIM>& G)
   {
      auto P = compute_P(G);
      auto A = compute_tangent_moduli(G);
      
      return std::make_pair(P, A);
   }

};
