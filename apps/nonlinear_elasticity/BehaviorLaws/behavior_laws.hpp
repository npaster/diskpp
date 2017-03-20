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

#include "material_parameters.hpp"
#include "common/eigen.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

/*

Fichier pour gérer les lois de comportements

-faire une boite qui prends en entrée eps --> behaviorbox --> stress
-penser aux variables internes pour sauvegarder partie elastique
- utilise t on des tenseur pour le module tangent
*/

/* Material: St Venant - Kirchhoof
 * Energy :  W(E) = lambda / 2 *(tr E)**2 + mu * tr (E**2)
 * Stress :  S(E) = lmabda * tr(E) * I_2 + 2.0*mu*E
 * Module :  C(E) = lambda * ( I_2 \cirl I_2)  + 2.0 * mu * I_4
 */


template<typename scalar_type>
class StVenantKirchhoffLaw
{
   typedef static_tensor4<scalar_type, 2, 2, 2, 2>   tensor_2d;
   typedef static_tensor4<scalar_type, 3, 3, 3, 3>   tensor_3d;
   scalar_type m_mu;
   scalar_type m_lambda;

   tensor_2d  elasticity_tensor_2d;
   tensor_3d  elasticity_tensor_3d;

   void
   init_elasticity_tensor()
   {
      elasticity_tensor_2d.setZero();

      for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
      for (int l = 0; l < 2; l++) {
         if(i == j && k == l)
         elasticity_tensor_2d(i,j,k,l) += m_lambda;

         if(i == k && j == l)
         elasticity_tensor_2d(i,j,k,l) += m_mu;

         if(i == l && j == k)
         elasticity_tensor_2d(i,j,k,l) += m_mu;
      }

      elasticity_tensor_3d.setZero();

      for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++) {
         if(i == j && k == l)
         elasticity_tensor_2d(i,j,k,l) += m_lambda;

         if(i == k && j == l)
         elasticity_tensor_2d(i,j,k,l) += m_mu;

         if(i == l && j == k)
         elasticity_tensor_2d(i,j,k,l) += m_mu;
      }
   }

public:
   StVenantKirchhoffLaw()
   : m_mu(0.0), m_lambda(0.0)
   {
      elasticity_tensor_2d.setZero();
      elasticity_tensor_3d.setZero();
   }

   StVenantKirchhoffLaw(const scalar_type mu, const scalar_type lambda)
   : m_mu(mu), m_lambda(lambda)
   {
      init_elasticity_tensor();
   }

   StVenantKirchhoffLaw(const Material_parameters<scalar_type>& material)
   : m_mu(material.giveMu()), m_lambda(material.giveLambda())
   {
      init_elasticity_tensor();
   }

   void
   setMu(const scalar_type mu)
   {
      m_mu = mu;
      init_elasticity_tensor();
   }

   void
   setLambda(const scalar_type lambda)
   {
      m_lambda = lambda;
      init_elasticity_tensor();
   }

   scalar_type
   giveMu() const {return m_mu;}

   scalar_type
   giveLambda() const {return m_lambda;}

   static_matrix<scalar_type,2,2>
   compute_PK2(const static_matrix<scalar_type,2,2>& CauchyGreenDroit)
   {
      static_matrix<scalar_type,2,2> Id = static_matrix<scalar_type,2,2>::Identity();
      static_matrix<scalar_type,2,2> EGl = 0.5 * (CauchyGreenDroit - Id);
      static_matrix<scalar_type,2,2> PK2 = m_lambda * EGl.trace() * Id + 2.0 * m_mu * EGl;

      return PK2;
   }

   static_matrix<scalar_type,3,3>
   compute_PK2(const static_matrix<scalar_type,3,3>& CauchyGreenDroit)
   {
      static_matrix<scalar_type,3,3> Id = static_matrix<scalar_type,3,3>::Identity();
      static_matrix<scalar_type,3,3> EGl = 0.5 * (CauchyGreenDroit - Id);
      static_matrix<scalar_type,3,3> PK2 = m_lambda * EGl.trace() * Id + 2.0 * m_mu * EGl;

      return PK2;
   }

   tensor_2d
   compute_tangent_moduli(const static_matrix<scalar_type,2,2>& CauchyGreenDroit)
   {
      return elasticity_tensor_2d;
   }

   tensor_3d
   compute_tangent_moduli(const static_matrix<scalar_type,3,3>& CauchyGreenDroit)
   {
      return elasticity_tensor_3d;
   }

};




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
class NeoNookeanLaw
{
   typedef static_tensor4<scalar_type, 2, 2, 2, 2>   tensor_2d;
   typedef static_tensor4<scalar_type, 3, 3, 3, 3>   tensor_3d;
   scalar_type m_mu;
   scalar_type m_lambda;
   size_t m_type;

   const size_t maxtype = 2;

   scalar_type
   computeU(scalar_type J)
   {
      if(m_type == 1)
         return ln(J);
      else if(m_type == 2)
         return sqrt(J);
      else
         assert(false);
   }

   scalar_type
   computeUprime(scalar_type J)
   {
      if(m_type == 1)
         return ln(J);
      else if(m_type == 2)
         return sqrt(J);
      else
         assert(false);
   }

public:
   NeoNookeanLaw()
   : m_mu(0.0), m_lambda(0.0), m_type(1)
   {
      if(m_type <= 0 || m_type > maxtype)
      {
         std::cout << "Unknown option for NeoNookean material" << '\n';
         std::cout << "We use U(J) = ln(J)" << '\n';
      }
   }

   NeoNookeanLaw(const scalar_type mu, const scalar_type lambda, const size_t type)
   : m_mu(mu), m_lambda(lambda), m_type(type)
   {
      if(m_type <= 0 || m_type > maxtype)
      {
         std::cout << "Unknown option for NeoNookean material" << '\n';
         std::cout << "We use U(J) = ln(J)" << '\n';
      }
   }

   NeoNookeanLaw(const Material_parameters<scalar_type>& material, const size_t type)
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

   scalar_type
   giveMu() const {return m_mu;}

   scalar_type
   giveLambda() const {return m_lambda;}

   static_matrix<scalar_type,2,2>
   compute_PK2(const static_matrix<scalar_type,2,2>& CauchyGreenDroit)
   {
      static_matrix<scalar_type,2,2> Id = static_matrix<scalar_type,2,2>::Identity();
      static_matrix<scalar_type,2,2> invCtr = (CauchyGreenDroit.inverse()).transpose();
      scalar_type J = sqrt(CauchyGreenDroit.det());
      scalar_type UJ = computeU(J);
      scalar_type UJp = computeUprime(J);

      return m_mu * (Id - invCtr) + m_lambda * J * UJ * UJp * invCtr;
   }

   static_matrix<scalar_type,3,3>
   compute_PK2(const static_matrix<scalar_type,3,3>& CauchyGreenDroit)
   {
      static_matrix<scalar_type,2,2> Id = static_matrix<scalar_type,2,2>::Identity();
      static_matrix<scalar_type,2,2> invCtr = (CauchyGreenDroit.inverse()).transpose();
      scalar_type J = sqrt(CauchyGreenDroit.det());
      scalar_type UJ = computeU(J);
      scalar_type UJp = computeUprime(J);

      return m_mu * (Id - invCtr) + m_lambda * J * UJ * UJp * invCtr;
   }

   tensor_2d
   compute_tangent_moduli(const static_matrix<scalar_type,2,2>& CauchyGreenDroit)
   {
      tensor_2d tangent_moduli = tensor_2d::Zero();
      return tangent_moduli;
   }

   tensor_3d
   compute_tangent_moduli(const static_matrix<scalar_type,3,3>& CauchyGreenDroit)
   {
      tensor_3d tangent_moduli = tensor_3d::Zero();
      return tangent_moduli;
   }

};
