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

#include "BehaviorLaws/material_parameters.hpp"
#include "common/eigen.hpp"
#include "BehaviorLaws/maths_tensor.hpp"
#include "BehaviorLaws/maths_deformation_tensors.hpp"
#include "BehaviorLaws/maths_jacobian_derivate.hpp"


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
   typedef scalar_type   tensor_1d;
   typedef static_tensor<scalar_type, 2>   tensor_2d;
   typedef static_tensor<scalar_type, 3>   tensor_3d;
   scalar_type m_mu;
   scalar_type m_lambda;

   tensor_1d  elasticity_tensor_1d;
   tensor_2d  elasticity_tensor_2d;
   tensor_3d  elasticity_tensor_3d;

   void
   init_elasticity_tensor()
   {
      elasticity_tensor_1d = m_lambda + 2.0 * m_mu;
      elasticity_tensor_2d = m_lambda * compute_IxI<scalar_type, 2>() + 2.0*m_mu * compute_IdentityTensor<scalar_type,2>();
      elasticity_tensor_3d = m_lambda * compute_IxI<scalar_type, 3>() + 2.0*m_mu * compute_IdentityTensor<scalar_type,3>();

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


   template<int DIM>
   static_matrix<scalar_type, DIM, DIM>
   compute_PK2(const static_matrix<scalar_type, DIM, DIM>& CauchyGreenDroit)
   {
      static_matrix<scalar_type, DIM, DIM> Id = static_matrix<scalar_type, DIM , DIM>::Identity();
      static_matrix<scalar_type, DIM, DIM> EGl = compute_GreenLagrangeTensor(CauchyGreenDroit);

      return m_lambda * EGl.trace() * Id + 2.0 * m_mu * EGl;
   }
   
   //1D case
   scalar_type
   compute_PK2(const scalar_type& CauchyGreenDroit)
   {
      scalar_type EGl = compute_GreenLagrangeTensor(CauchyGreenDroit);
      
      return m_lambda * EGl  + 2.0 * m_mu * EGl;
   }

   tensor_1d
   compute_tangent_moduli(const static_matrix<scalar_type,1,1>& CauchyGreenDroit)
   {
      return elasticity_tensor_1d;
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

   std::pair<scalar_type, scalar_type>
   compute_whole(const scalar_type& CauchyGreenDroit)
   {
      scalar_type PK2 = compute_PK2(CauchyGreenDroit);
      
      return std::make_pair(PK2, elasticity_tensor_1d);
   }

   std::pair<static_matrix<scalar_type, 2, 2>, static_tensor<scalar_type, 2> >
   compute_whole(const static_matrix<scalar_type, 2, 2>& CauchyGreenDroit)
   {
      static_matrix<scalar_type, 2, 2> PK2 = compute_PK2(CauchyGreenDroit);

      return std::make_pair(PK2, elasticity_tensor_2d);
   }

   std::pair<static_matrix<scalar_type, 3, 3>, static_tensor<scalar_type, 3> >
   compute_whole(const static_matrix<scalar_type, 3, 3>& CauchyGreenDroit)
   {
      static_matrix<scalar_type, 3, 3> PK2 = compute_PK2(CauchyGreenDroit);

      return std::make_pair(PK2, elasticity_tensor_3d);
   }

};
