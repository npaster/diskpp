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
#include "bases/bases_utils.hpp"

// the fourth order tensor is stored like a Matrix
// | A1111  A1112  A1211  A1212 |
// | A1121  A1122  A1221  A1222 |
// | A2111  A2112  A2211  A2212 |
// | A2112  A2122  A2221  A2222 |


// Convert Matrix in Coulum vector

template<typename T, int DIM>
void
converttovector(const static_matrix<T, DIM, DIM>& mat)
{
   static_assert( (DIM == 2 || DIM == 3), "Can not compute conversion for this dimension");
}



template<typename T>
static_vector<T, 4>
converttovector(const static_matrix<T, 2, 2>& mat)
{
   return static_vector<T,4>{mat(0,0), mat(0,1), mat(1,0), mat(1,1)};
}

template<typename T>
static_vector<T, 9>
converttovector(const static_matrix<T, 3, 3>& mat)
{
   static_vector<T,9> ret;

   ret(0) = mat(0,0);
   ret(1) = mat(0,1);
   ret(2) = mat(0,2);

   ret(3) = mat(1,0);
   ret(4) = mat(1,1);
   ret(4) = mat(1,2);

   ret(6) = mat(2,0);
   ret(7) = mat(2,1);
   ret(8) = mat(2,2);

   return ret;
}


// Convert vector in matrix

template<typename T, int DIM>
void
converttomatrix(const static_vector<T, DIM>& vec)
{
   static_assert( (DIM == 4 || DIM == 9), "Can not compute conversion for this dimension");
}



template<typename T>
static_matrix<T, 2, 2>
converttomatrix(const static_vector<T, 4>& vec)
{
   static_matrix<T, 2, 2> mat;

   mat(0,0) = vec(0);
   mat(0,1) = vec(1);

   mat(1,0) = vec(2);
   mat(1,1) = vec(1);

   return mat;
}

template<typename T>
static_matrix<T, 3, 3>
converttomatrix(const static_vector<T, 9>& vec)
{
   static_matrix<T, 2, 2> mat;

   mat(0,0) = vec(0);
   mat(0,1) = vec(1);
   mat(0,2) = vec(2);

   mat(1,0) = vec(3);
   mat(1,1) = vec(4);
   mat(1,2) = vec(5);

   mat(2,0) = vec(6);
   mat(2,1) = vec(7);
   mat(2,2) = vec(8);

   return mat;
}



// Product Tensor - Matrix

template<typename T, int  DIM>
static_matrix<T, DIM, DIM>
tm_prod(const static_tensor<T, DIM>& tens, const static_matrix<T, DIM, DIM>& mat )
{
   static_matrix<T, DIM, DIM> ret;

   for (size_t i = 0; i < DIM; i++)
      for (size_t j = 0; j < DIM; j++)
         ret(i,j) = mm_prod(tens.block(i*DIM, j*DIM, DIM, DIM), mat);

   return ret;
}


// Compute Kronecker product

template<typename T, int M, int N, int P, int Q>
void
computeKroneckerProduct(const static_matrix<T, M, N>& A, const static_matrix<T, P, Q>& B )
{
   static_assert( (M == N && N == P && P == Q),  "Kronecker product : Not yet develloped");
}

template<typename T, int DIM>
static_tensor<T, DIM>
computeKroneckerProduct(const static_matrix<T, DIM, DIM>& A, const static_matrix<T, DIM, DIM>& B )
{
   static_tensor<T, DIM> ret;

   for (size_t i = 0; i < DIM; i++)
      for (size_t j = 0; j < DIM; j++)
         ret.block(i*DIM, j*DIM, DIM, DIM) = A(i,j) * B;

   return ret;
}


template<typename T, int  DIM>
static_tensor<T, DIM>
compute_IdentityTensor()
{
   static_tensor<T, DIM> ret = static_tensor<T,DIM>::Zero();
   T one = T{1};

   if(DIM == 2){
      ret(0, 0) = one;  // I1111
      ret(0, 3) = one;  // I1212
      ret(3, 0) = one;  // I2121
      ret(3, 3) = one;  // I2222
   }
   else if( DIM == 3){
      ret(0, 0) = one;  // I1111
      ret(0, 4) = one;  // I1212
      ret(0, 8) = one;  // I1313
      ret(4, 0) = one;  // I2121
      ret(4, 4) = one;  // I2222
      ret(4, 8) = one;  // I2323
      ret(8, 0) = one;  // I3131
      ret(8, 4) = one;  // I3232
      ret(8, 8) = one;  // I3333
   }
   else
      static_assert((DIM == 2 || DIM == 3), "Wrong dimension only 2 and 3");

   return ret;
}


template<typename T, int DIM>
static_tensor<T, DIM>
compute_IxI()
{
   return static_tensor<T,DIM>::Identity();
}
