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

#pragma once

#include "common/eigen.hpp"

// Function to convert Stress to an other
// F is the deformation gradient

namespace disk {

// PK2 = F^{-1} * PK1 

template<typename T, int  DIM>
static_matrix<T, DIM, DIM>
convertPK1toPK2(const static_matrix<T, DIM, DIM>& PK1, const static_matrix<T, DIM, DIM>& F)
{
   return (F.inverse()) * PK1;
}

// PK1 = F * PK2
template<typename T, int  DIM>
static_matrix<T, DIM, DIM>
convertPK2toPK1(const static_matrix<T, DIM, DIM>& PK2, const static_matrix<T, DIM, DIM>& F)
{
   return F * PK2;
}

// Cauchy = J^{-1} * P * F^{T}
template<typename T, int  DIM>
static_matrix<T, DIM, DIM>
convertPK1toCauchy(const static_matrix<T, DIM, DIM>& PK1, const static_matrix<T, DIM, DIM>& F)
{
   return (PK1 * F.transpose())/ F.determinant();
}

// PK1 = J * Cauchy * F^{-T}
template<typename T, int  DIM>
static_matrix<T, DIM, DIM>
convertCauchytoPK1(const static_matrix<T, DIM, DIM>& Cauchy, const static_matrix<T, DIM, DIM>& F)
{
   const auto invF = F.inverse();
   return F.determinant() * Cauchy * invF.transpose();
}

// PK2 = J  * F^{-1} * Cauchy * F^{-T}
template<typename T, int  DIM>
static_matrix<T, DIM, DIM>
convertCauchytoPK2(const static_matrix<T, DIM, DIM>& Cauchy, const static_matrix<T, DIM, DIM>& F)
{
   const auto invF = F.inverse();
   return F.determinant() * invF * Cauchy * invF.transpose();
}


// Cauchy = J^-1 * F * Cauchy * F^{T}
template<typename T, int  DIM>
static_matrix<T, DIM, DIM>
convertPK2toCauchy(const static_matrix<T, DIM, DIM>& PK2, const static_matrix<T, DIM, DIM>& F)
{
   return  F * PK2 * F.transpose()/F.determinant();
}

}//end namespace