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
    scalar_type m_mu;
    scalar_type m_lambda;


 public:
     StVenantKirchhoffLaw()
     : m_mu(0.0), m_lambda(0.0) {}

     StVenantKirchhoffLaw(const scalar_type mu, const scalar_type lambda)
     : m_mu(mu), m_lambda(lambda) {}

     StVenantKirchhoffLaw(const Material_parameters<scalar_type>& material)
     : m_mu(material.giveMu()), m_lambda(material.giveLambda()) {}

    void
    setMu(const scalar_type mu) { m_mu = mu;}

    void
    setLambda(const scalar_type lambda) { m_lambda = lambda;}

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

};
