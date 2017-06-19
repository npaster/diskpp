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

// Compute F = G + I

template<typename T>
T
compute_FTensor(const T& Gradient)
{
   return Gradient + T{1};
}

template<typename T, int  DIM>
static_matrix<T, DIM, DIM>
compute_FTensor(const static_matrix<T, DIM, DIM>& Gradient)
{
   return Gradient + static_matrix<T,DIM, DIM>::Identity();
}

// Compute C = F^T * F

template<typename T>
T
compute_CauchyGreenRightTensor(const T& FTensor)
{
   return FTensor * FTensor;
}

template<typename T, int  DIM>
static_matrix<T, DIM, DIM>
compute_CauchyGreenRightTensor(const static_matrix<T, DIM, DIM>& FTensor)
{
   return FTensor.transpose() * FTensor;
}

//Compute E = 1/2 *( C - I)

template<typename T>
T
compute_GreenLagrangeTensor(const T& CauchyGreenTensor)
{
   return 0.5 * (CauchyGreenTensor - T{1});
}


template<typename T, int  DIM>
static_matrix<T, DIM, DIM>
compute_GreenLagrangeTensor(const static_matrix<T, DIM, DIM>& CauchyGreenTensor)
{
   return 0.5 * (CauchyGreenTensor - static_matrix<T,DIM, DIM>::Identity());
}

// Compute b = F * F^T

template<typename T, int  DIM>
static_matrix<T, DIM, DIM>
compute_CauchyGreenLeftTensor(const static_matrix<T, DIM, DIM>& FTensor)
{
   return FTensor * FTensor.transpose();
}


//convert C to A

template<typename T, int  DIM>
static_tensor<T, DIM>
convertto_elastic_moduli_2(const static_matrix<T, DIM, DIM>& F, const static_matrix<T, DIM, DIM>& PK2,
                         const static_tensor<T, DIM>& C)
{
   static_matrix<T, DIM, DIM> Id = static_matrix<T, DIM,DIM>::Identity();
   return computeKroneckerProduct( PK2, Id) + 2.0 * computeKroneckerProduct(tm_prod(C,F), F );
}
