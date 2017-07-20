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



/* Material: Neo-nookean
 * Energy :  W(F) = Wiso(F) + Wvol(F)
 *   - Wiso(F) =    mu / 2 *[tr(F^T * F) - d] - mu * ln(J)
 *   - Wvol(F) =  lambda/2 * U(J)**2
 * ** We set T1(J) = J * U(J) * U'(J) and T2(J) =  U(J) * J *( U''(J) * J + U'(J)) + ( J * U'(J))^2
 * Stress :  PK1(F) = Piso(F) + Pvol(F)
 *   - Piso(F) = mu * ( F - F^{-T})
 *   - Pvol(F) = lambda * T1(J) * F^{-T}
 * Module :  A(F) = PK1(F) = Aiso(F) + Avol(F)
 *   - Aiso(F) = mu * ( I4 + F^{-T} \time_inf F^{-1})
 *   - Avol(F) = -lambda * T1(J) * F^{-T} \time_inf F^{-1}
 *                  + lambda * T2(J) * F^{-T} \kronecker F^{-T}
 */

/* Laws:
 * 1- U(J) = ln(J)
 * 2- U(J) = J -1
 * 3- U(J) = log10(J)
 * 4- U(J) = 1 -1 /J
 * 5- U(J) = J^2 -1
 * 6- U(J) = sqrt( ( J^2 -1 - 2 *ln(J)) /2)
 * */


template<typename scalar_type>
class NeoHookeanLaw
{
   scalar_type m_mu;
   scalar_type m_lambda;
   size_t m_type;

   const size_t maxtype = 6;

   scalar_type
   compute_T1(scalar_type J)
   {
      if(m_type == 1)
         return log(J);
      else if(m_type == 2)
         return J * ( J - 1.0);
      else if(m_type == 3)
         return log(J)/(log(10)*log(10));
      else if(m_type == 4)
         return (J - 1.0)/(J*J);
      else if(m_type == 5){
         scalar_type J2 = J *J;
         return 2*J2*(J2 -1.0);
      }
      else if(m_type == 6)
         return (J*J -1.0)/2.0;
      else
         throw std::invalid_argument("NeoHookeanLaw: m_type have to be <= 6");
   }

   scalar_type
   compute_T2(scalar_type J)
   {
      if(m_type == 1)
         return 1.0;
      else if(m_type == 2)
         return J * ( 2.0 * J - 1.0);
      else if(m_type == 3)
         return 1.0/(log(10)*log(10));
      else if(m_type == 4)
         return (2.0 - J)/(J*J);
      else if(m_type == 5){
         scalar_type J2 = J *J;
         return J2*(8.0*J2 -4.0);
      }
      else if(m_type == 6)
         return J*J;
      else
         throw std::invalid_argument("NeoHookeanLaw: m_type have to be <= 6");
   }


public:
   NeoHookeanLaw()
   : m_mu(0.0), m_lambda(0.0), m_type(1)
   {
      if(m_type <= 0 || m_type > maxtype)
      {
         std::cout << "Unknown option for NeoNookean material" << '\n';
         std::cout << "We use U(J) = ln(J)" << '\n';
         m_type = 1;
      }
   }

   NeoHookeanLaw(const scalar_type mu, const scalar_type lambda, const size_t type)
   : m_mu(mu), m_lambda(lambda), m_type(type)
   {
      if(m_type <= 0 || m_type > maxtype)
      {
         std::cout << "Unknown option for NeoNookean material" << '\n';
         std::cout << "We use U(J) = ln(J)" << '\n';
         m_type = 1;
      }
   }

   void
   setMu(const scalar_type mu)
   {
      m_mu = mu;
   }

   void
   setLambda(const scalar_type lambda)
   {
      m_lambda = lambda;
   }

   void
   setType(const size_t type)
   {
      m_type = type;
      if(m_type <= 0 || m_type > maxtype)
      {
         std::cout << "Unknown option for NeoNookean material" << '\n';
         std::cout << "We use U(J) = ln(J)" << '\n';
         m_type = 1;
      }
   }

   scalar_type
   giveMu() const {return m_mu;}

   scalar_type
   giveLambda() const {return m_lambda;}


   template< int DIM>
   static_matrix<scalar_type, DIM, DIM>
   compute_PK1(const static_matrix<scalar_type, DIM, DIM>& F)
   {
      static_matrix<scalar_type, DIM, DIM> invF = F.inverse();
      scalar_type J = F.determinant();
      scalar_type T1 = compute_T1(J);

      if(J <=0.0)
         throw std::invalid_argument("J <= 0");

      return  m_mu * F + ( m_lambda * T1 - m_mu) * invF.transpose();
   }


   template<int DIM>
   static_tensor<scalar_type, DIM>
   compute_tangent_moduli_A(const static_matrix<scalar_type, DIM, DIM>& F)
   {
      static_matrix<scalar_type, DIM, DIM> invF = F.inverse();
      static_matrix<scalar_type, DIM, DIM> invFt = invF.transpose();
      scalar_type J = F.determinant();
      scalar_type T1 = compute_T1(J);
      scalar_type T2 = compute_T2(J);
      
      static_tensor<scalar_type, DIM> I4 = compute_IdentityTensor<scalar_type,DIM>();
      static_tensor<scalar_type, DIM> invFt_invF = computeProductInf(invFt, invF);
      static_tensor<scalar_type, DIM> invFt_invFt = computeKroneckerProduct(invFt, invFt);

      if(J <=0.0)
         throw std::invalid_argument("J <= 0");

      return m_mu * (I4 + invFt_invF) + m_lambda *( T2 * invFt_invFt - T1 * invFt_invF);
   }


   template<int DIM>
   std::pair<static_matrix<scalar_type, DIM, DIM>, static_tensor<scalar_type, DIM> >
   compute_whole_PK1(const static_matrix<scalar_type, DIM, DIM>& F)
   {
      static_matrix<scalar_type, DIM, DIM> invF = F.inverse();
      static_matrix<scalar_type, DIM, DIM> invFt = invF.transpose();

      scalar_type J = F.determinant();
      scalar_type T1 = compute_T1(J);
      scalar_type T2 = compute_T2(J);


      if(J <=0.0){
         std::cout << "J = " << J << std::endl;
         throw std::invalid_argument("J <= 0");
      }

      static_tensor<scalar_type, DIM> I4 = compute_IdentityTensor<scalar_type,DIM>();
      static_tensor<scalar_type, DIM> invFt_invF = computeProductInf(invFt, invF);
      static_tensor<scalar_type, DIM> invFt_invFt = computeKroneckerProduct(invFt, invFt);

      auto PK1 = m_mu * F + ( m_lambda * T1 - m_mu) * invFt;

      auto A = m_mu * (I4 + invFt_invF) + m_lambda *( T2 * invFt_invFt - T1 * invFt_invF);

      return std::make_pair(PK1, A);
   }

};
