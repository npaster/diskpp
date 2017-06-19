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

#include "BehaviorLaws/material_parameters.hpp"
#include "common/eigen.hpp"
#include "BehaviorLaws/maths_tensor.hpp"
#include "BehaviorLaws/maths_jacobian_derivate.hpp"
#include "BehaviorLaws/maths_deformation_tensors.hpp"
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
 * Energy :  W(C) = Wiso(C) + Wvol(C)
 *   - Wiso(C) =    mu / 2 *[tr(C) - d] - mu * ln(J)
 *   - Wvol(C) =  lambda/2 * U(J)**2
 * Stress :  S(C) = Siso(C) + Svol(C)
 *   - Siso(C) = mu * ( Id - C^{-T})
 *   - Svol(C) = lambda * J * U(J) * U'(J)
 * Module :  C(E) = lambda * ( I_2 \cirl I_2)  + 2.0 * mu * I_4
 */


template<typename scalar_type>
class NeoHookeanLaw
{
   scalar_type m_mu;
   scalar_type m_lambda;
   size_t m_type;

   const size_t maxtype = 3;

   scalar_type
   computeU(scalar_type J)
   {
      if(m_type == 1)
         return J - scalar_type{1};
      else if(m_type == 2)
         return log(J);
      else if(m_type == 3)
         return scalar_type{1} - scalar_type{1}/J;
      else
         assert(false);
   }

   scalar_type
   computeUprime(scalar_type J)
   {
      if(m_type == 1)
         return scalar_type{1};
      else if(m_type == 2)
         return scalar_type{1}/J;
      else if(m_type == 3)
         return -scalar_type{1}/(J*J);
      else
         assert(false);
   }

   scalar_type
   computeUsecond(scalar_type J)
   {
      if(m_type == 1)
         return scalar_type{0};
      else if(m_type == 2)
         return -scalar_type{1}/(J*J);
      else if(m_type == 3)
         return scalar_type{2}/(J*J*J);
      else
         assert(false);
   }

public:
   NeoHookeanLaw()
   : m_mu(0.0), m_lambda(0.0), m_type(1)
   {
      if(m_type <= 0 || m_type > maxtype)
      {
         std::cout << "Unknown option for NeoNookean material" << '\n';
         std::cout << "We use U(J) = ln(J)" << '\n';
      }
   }

   NeoHookeanLaw(const scalar_type mu, const scalar_type lambda, const size_t type)
   : m_mu(mu), m_lambda(lambda), m_type(type)
   {
      if(m_type <= 0 || m_type > maxtype)
      {
         std::cout << "Unknown option for NeoNookean material" << '\n';
         std::cout << "We use U(J) = ln(J)" << '\n';
      }
   }

   NeoHookeanLaw(const Material_parameters<scalar_type>& material, const size_t type)
   : m_mu(material.giveMu()), m_lambda(material.giveLambda()), m_type(type)
   {
      if(m_type <= 0 || m_type > maxtype)
      {
         std::cout << "Unknown option for NeoNookean material" << '\n';
         std::cout << "We use U(J) = ln(J)" << '\n';
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
   }

   scalar_type
   giveMu() const {return m_mu;}

   scalar_type
   giveLambda() const {return m_lambda;}

   template< int DIM>
   static_matrix<scalar_type, DIM, DIM>
   compute_PK2(const static_matrix<scalar_type, DIM, DIM>& CauchyGreenDroit)
   {
      static_matrix<scalar_type, DIM, DIM> Id = static_matrix<scalar_type, DIM, DIM>::Identity();
      static_matrix<scalar_type, DIM, DIM> dJ_dC = computeJacobianFirstDerivate(CauchyGreenDroit);
      scalar_type J = sqrt(CauchyGreenDroit.determinant());
      scalar_type UJ = computeU(J);
      scalar_type UJp = computeUprime(J);

      auto Siso = m_mu * (Id - scalar_type{2.0} * dJ_dC / J);
      auto Svol = scalar_type{2.0} * m_lambda * UJ * UJp * dJ_dC;

      return  Siso + Svol;
   }



   template<int DIM>
   static_tensor<scalar_type, DIM>
   compute_tangent_moduli_C(const static_matrix<scalar_type, DIM, DIM>& CauchyGreenDroit)
   {
      static_matrix<scalar_type, DIM, DIM> dJ_dC = computeJacobianFirstDerivate(CauchyGreenDroit);
      static_tensor<scalar_type, DIM> d2J_dC2 = computeJacobianSecondDerivate(CauchyGreenDroit);
      static_tensor<scalar_type, DIM> square_dJ_dC = computeKroneckerProduct(dJ_dC, dJ_dC);
      scalar_type J = sqrt(CauchyGreenDroit.determinant());
      scalar_type UJ = computeU(J);
      scalar_type UJp = computeUprime(J);
      scalar_type UJs = computeUsecond(J);

      auto Kiso = - scalar_type{4} * m_mu * ( - square_dJ_dC / (J*J) + d2J_dC2 / J);
      auto Kvol = scalar_type{4} * m_lambda * ((UJp * UJp + UJ * UJs) * square_dJ_dC + UJ * UJp * d2J_dC2);

      return Kiso + Kvol;
   }


   template<int DIM>
   std::pair<static_matrix<scalar_type, DIM, DIM>, static_tensor<scalar_type, DIM> >
   compute_whole_PK2(const static_matrix<scalar_type, DIM, DIM>& CauchyGreenDroit)
   {
      static_matrix<scalar_type, DIM, DIM> Id = static_matrix<scalar_type, DIM, DIM>::Identity();
      static_matrix<scalar_type, DIM, DIM> dJ_dC = computeJacobianFirstDerivate(CauchyGreenDroit);
      static_tensor<scalar_type, DIM> d2J_dC2 = computeJacobianSecondDerivate(CauchyGreenDroit);
      static_tensor<scalar_type, DIM> square_dJ_dC = computeKroneckerProduct(dJ_dC, dJ_dC);

      scalar_type J = sqrt(CauchyGreenDroit.determinant());
      scalar_type UJ = computeU(J);
      scalar_type UJp = computeUprime(J);
      scalar_type UJs = computeUsecond(J);

      auto Siso = m_mu * (Id - scalar_type{2.0} * dJ_dC / J);
      auto Svol = scalar_type{2.0} * m_lambda * UJ * UJp * dJ_dC;

      auto Kiso = - scalar_type{4} * m_mu * ( - square_dJ_dC / (J*J) + d2J_dC2 / J);
      auto Kvol = scalar_type{4} * m_lambda * ((UJp * UJp + UJ * UJs) * square_dJ_dC + UJ * UJp * d2J_dC2);

      return std::make_pair(Siso + Svol, Kiso + Kvol);
   }
   
   
   template< int DIM>
   static_matrix<scalar_type, DIM, DIM>
   compute_PK1(const static_matrix<scalar_type, DIM, DIM>& F)
   {
      static_matrix<scalar_type, DIM, DIM> invF = F.inverse();
      scalar_type J = F.determinant();
      
      return  m_mu * F + ( m_lambda * log(J) - m_mu) * invF.transpose();
   }
   
   
   
   template<int DIM>
   static_tensor<scalar_type, DIM>
   compute_tangent_moduli_A(const static_matrix<scalar_type, DIM, DIM>& F)
   {
      static_matrix<scalar_type, DIM, DIM> invF = F.inverse();
      static_matrix<scalar_type, DIM, DIM> invFt = invF.transpose();
      scalar_type J = F.determinant();
      
      static_tensor<scalar_type, DIM> I4 = compute_IdentityTensor<scalar_type,DIM>();
      
      return m_lambda * computeKroneckerProduct(invFt, invFt) + m_mu * I4 + (m_mu - m_lambda * log(J)) * computeProductInf(invFt, invF);
   }
   
   
   template<int DIM>
   std::pair<static_matrix<scalar_type, DIM, DIM>, static_tensor<scalar_type, DIM> >
   compute_whole_PK1(const static_matrix<scalar_type, DIM, DIM>& F)
   {
      static_matrix<scalar_type, DIM, DIM> invF = F.inverse();
      static_matrix<scalar_type, DIM, DIM> invFt = invF.transpose();
      scalar_type J = F.determinant();
      
      if(J <=0.0) 
         std::cout << "J <= 0 : " << J << std::endl;
      
      static_tensor<scalar_type, DIM> I4 = compute_IdentityTensor<scalar_type,DIM>();
      
      auto PK1 = m_mu * F + ( m_lambda * log(J) - m_mu) * invF.transpose();
      
      auto A = m_lambda * computeKroneckerProduct(invFt, invFt) + m_mu * I4 + (m_mu - m_lambda * log(J)) * computeProductInf(invFt, invF);
      
      return std::make_pair(PK1, A);
   }

};
