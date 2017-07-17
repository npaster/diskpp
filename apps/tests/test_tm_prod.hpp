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

#include <iostream>

#include "common/eigen.hpp"
#include "BehaviorLaws/maths_tensor.hpp"

#pragma once


template< typename T>
void
test_tm_prod(const T prec)
{
   static_matrix<T, 3, 3> G;
   G(0,0) = T{1.0}; G(0,1) = T{2.0}; G(0,2) = T{3.0};
   G(1,0) = T{4.0}; G(1,1) = T{6.0}; G(1,2) = T{8.0};
   G(2,0) = T{5.0}; G(2,1) = T{7.0}; G(2,2) = T{9.0};

   static_matrix<T, 3, 3> R;
   R(0,0) = T{38.0}; R(0,1) = T{9.0}; R(0,2) = T{67.0};
   R(1,0) = T{48.0}; R(1,1) = T{52.0}; R(1,2) = T{63.0};
   R(2,0) = T{13.0}; R(2,1) = T{0}; R(2,2) = T{28.0};

   static_tensor<T, 3> A = static_tensor<T, 3>::Zero();
   A(0,0) = T{1.0}; A(0,3) = T{2.0}; A(0,8) = T{6.0};
   A(1,2) = T{2.0}; A(1,7) = T{4.0};
   A(2,1) = T{3.0}; A(2,4) = T{1.0}; A(2,6) = T{5.0};
   A(3,1) = T{3.0}; A(3,3) =-T{1.0}; A(3,7) = T{3.0};
   A(4,2) = T{4.0}; A(4,3) = T{2.0}; A(4,6) = T{3.0}; A(4,8) = T{3.0};
   A(5,0) = T{2.0}; A(5,5) = T{5.0}; A(5,7) = T{3.0};
   A(6,0) = T{1.0};
   A(7,1) = T{2.0}; A(7,6) = T{5.0};A(7,8) = T{1.0};


   T error_R = (R - tm_prod(A,G)).norm();

   if(error_R > prec){
      std::cout << "tm_prod, error P: " << error_R << std::endl;
      std::cout << "P analytique: " << R << std::endl;
      std::cout << "P calculé: " << tm_prod(A,G) << std::endl;
   }
}


template< typename T>
void
test_kronecker(const T prec)
{
   static_matrix<T, 3, 3> G;
   G(0,0) = T{1.0}; G(0,1) = T{2.0}; G(0,2) = T{3.0};
   G(1,0) = T{4.0}; G(1,1) = T{6.0}; G(1,2) = T{8.0};
   G(2,0) = T{5.0}; G(2,1) = T{7.0}; G(2,2) = T{9.0};

   static_matrix<T, 3, 3> R;
   R(0,0) = T{1.0}; R(0,1) = T{9.0}; R(0,2) = T{6.0};
   R(1,0) = T{8.0}; R(1,1) = T{2.0}; R(1,2) = T{3.0};
   R(2,0) = T{3.0}; R(2,1) = T{0};   R(2,2) = T{2.0};

   static_tensor<T, 3> A = static_tensor<T, 3>::Zero();
   //A11..
   A(0,0) = T{1.0}; A(0,1) = T{9.0}; A(0,2) = T{6.0};
   A(1,0) = T{8.0}; A(1,1) = T{2.0}; A(1,2) = T{3.0};
   A(2,0) = T{3.0}; A(2,1) = T{0.0}; A(2,2) = T{2.0};

   //A12..
   A(0,3) = T{2.0}; A(0,4) = T{18.0}; A(0,5) = T{12.0};
   A(1,3) = T{16.0}; A(1,4) = T{4.0}; A(1,5) = T{6.0};
   A(2,3) = T{6.0}; A(2,4) = T{0.0}; A(2,5) = T{4.0};

   //A13..
   A(0,6) = T{3.0}; A(0,7) = T{27.0}; A(0,8) = T{18.0};
   A(1,6) = T{24.0}; A(1,7) = T{6.0}; A(1,8) = T{9.0};
   A(2,6) = T{9.0}; A(2,7) = T{0.0}; A(2,8) = T{6.0};

   //A21..
   A(3,0) = T{4.0}; A(3,1) = T{36.0}; A(3,2) = T{24.0};
   A(4,0) = T{32.0}; A(4,1) = T{8.0}; A(4,2) = T{12.0};
   A(5,0) = T{12.0}; A(5,1) = T{0.0}; A(5,2) = T{8.0};

   //A22..
   A(3,3) = T{6.0}; A(3,4) = T{54.0}; A(3,5) = T{36.0};
   A(4,3) = T{48.0}; A(4,4) = T{12.0}; A(4,5) = T{18.0};
   A(5,3) = T{18.0}; A(5,4) = T{0.0}; A(5,5) = T{12.0};

   //A23..
   A(3,6) = T{8.0}; A(3,7) = T{72.0}; A(3,8) = T{48.0};
   A(4,6) = T{64.0}; A(4,7) = T{16.0}; A(4,8) = T{24.0};
   A(5,6) = T{24.0}; A(5,7) = T{0.0}; A(5,8) = T{16.0};

   //A31..
   A(6,0) = T{5.0}; A(6,1) = T{45.0}; A(6,2) = T{30.0};
   A(7,0) = T{40.0}; A(7,1) = T{10.0}; A(7,2) = T{15.0};
   A(8,0) = T{15.0}; A(8,1) = T{0.0}; A(8,2) = T{10.0};

   //A32..
   A(6,3) = T{7.0}; A(6,4) = T{63.0}; A(6,5) = T{42.0};
   A(7,3) = T{56.0}; A(7,4) = T{14.0}; A(7,5) = T{21.0};
   A(8,3) = T{21.0}; A(8,4) = T{0.0}; A(8,5) = T{14.0};

   //A33..
   A(6,6) = T{9.0}; A(6,7) = T{81.0}; A(6,8) = T{54.0};
   A(7,6) = T{72.0}; A(7,7) = T{18.0}; A(7,8) = T{27.0};
   A(8,6) = T{27.0}; A(8,7) = T{0.0}; A(8,8) = T{18.0};


   auto K = computeKroneckerProduct(G, R);
   T error_R = (A -  K).norm();

   if(error_R > prec){
      std::cout << "kronecker, error P: " << error_R << std::endl;
      std::cout << "P analytique: " << A << std::endl;
      std::cout << "P calculé: " << K << std::endl;
   }

}


template< typename T>
void
test_prodsup(const T prec)
{
   static_matrix<T, 3, 3> G;
   G(0,0) = T{1.0}; G(0,1) = T{2.0}; G(0,2) = T{3.0};
   G(1,0) = T{4.0}; G(1,1) = T{6.0}; G(1,2) = T{8.0};
   G(2,0) = T{5.0}; G(2,1) = T{7.0}; G(2,2) = T{9.0};

   static_matrix<T, 3, 3> R;
   R(0,0) = T{1.0}; R(0,1) = T{9.0}; R(0,2) = T{6.0};
   R(1,0) = T{8.0}; R(1,1) = T{2.0}; R(1,2) = T{3.0};
   R(2,0) = T{3.0}; R(2,1) = T{0};   R(2,2) = T{2.0};

   static_tensor<T, 3> A = static_tensor<T, 3>::Zero();
   //A11..
   A(0,0) = T{1.0}; A(0,1) = T{9.0}; A(0,2) = T{6.0};
   A(1,0) = T{2.0}; A(1,1) = T{18.0}; A(1,2) = T{12.0};
   A(2,0) = T{3.0}; A(2,1) = T{27.0}; A(2,2) = T{18.0};

   //A12..
   A(0,3) = T{8.0}; A(0,4) = T{2.0}; A(0,5) = T{3.0};
   A(1,3) = T{16.0}; A(1,4) = T{4.0}; A(1,5) = T{6.0};
   A(2,3) = T{24.0}; A(2,4) = T{6.0}; A(2,5) = T{9.0};

   //A13..
   A(0,6) = T{3.0}; A(0,7) = T{0.0}; A(0,8) = T{2.0};
   A(1,6) = T{6.0}; A(1,7) = T{0.0}; A(1,8) = T{4.0};
   A(2,6) = T{9.0}; A(2,7) = T{0.0}; A(2,8) = T{6.0};

   //A21..
   A(3,0) = T{4.0}; A(3,1) = T{36.0}; A(3,2) = T{24.0};
   A(4,0) = T{6.0}; A(4,1) = T{54.0}; A(4,2) = T{36.0};
   A(5,0) = T{8.0}; A(5,1) = T{72.0}; A(5,2) = T{48.0};

   //A22..
   A(3,3) = T{32.0}; A(3,4) = T{8.0}; A(3,5) = T{12.0};
   A(4,3) = T{48.0}; A(4,4) = T{12.0}; A(4,5) = T{18.0};
   A(5,3) = T{64.0}; A(5,4) = T{16.0}; A(5,5) = T{24.0};

   //A23..
   A(3,6) = T{12.0}; A(3,7) = T{0.0}; A(3,8) = T{8.0};
   A(4,6) = T{18.0}; A(4,7) = T{0.0}; A(4,8) = T{12.0};
   A(5,6) = T{24.0}; A(5,7) = T{0.0}; A(5,8) = T{16.0};

   //A31..
   A(6,0) = T{5.0}; A(6,1) = T{45.0}; A(6,2) = T{30.0};
   A(7,0) = T{7.0}; A(7,1) = T{63.0}; A(7,2) = T{42.0};
   A(8,0) = T{9.0}; A(8,1) = T{81.0}; A(8,2) = T{54.0};

   //A32..
   A(6,3) = T{40.0}; A(6,4) = T{10.0}; A(6,5) = T{15.0};
   A(7,3) = T{56.0}; A(7,4) = T{14.0}; A(7,5) = T{21.0};
   A(8,3) = T{72.0}; A(8,4) = T{18.0}; A(8,5) = T{27.0};

   //A33..
   A(6,6) = T{15.0}; A(6,7) = T{0.0}; A(6,8) = T{10.0};
   A(7,6) = T{21.0}; A(7,7) = T{0.0}; A(7,8) = T{14.0};
   A(8,6) = T{27.0}; A(8,7) = T{0.0}; A(8,8) = T{18.0};


   auto K = computeProductSup(G, R);
   T error_R = (A -  K).norm();

   if(error_R > prec){
      std::cout << "prodsup, error P: " << error_R << std::endl;
      std::cout << "P analytique: " << A << std::endl;
      std::cout << "P calculé: " << K << std::endl;
   }

}


template< typename T>
void
test_prodinf(const T prec)
{
   static_matrix<T, 3, 3> G;
   G(0,0) = T{1.0}; G(0,1) = T{2.0}; G(0,2) = T{3.0};
   G(1,0) = T{4.0}; G(1,1) = T{6.0}; G(1,2) = T{8.0};
   G(2,0) = T{5.0}; G(2,1) = T{7.0}; G(2,2) = T{9.0};

   static_matrix<T, 3, 3> R;
   R(0,0) = T{1.0}; R(0,1) = T{9.0}; R(0,2) = T{6.0};
   R(1,0) = T{8.0}; R(1,1) = T{2.0}; R(1,2) = T{3.0};
   R(2,0) = T{3.0}; R(2,1) = T{0};   R(2,2) = T{2.0};

   static_tensor<T, 3> A = static_tensor<T, 3>::Zero();
   //A11..
   A(0,0) = T{1.0}; A(0,1) = T{2.0}; A(0,2) = T{3.0};
   A(1,0) = T{9.0}; A(1,1) = T{18.0}; A(1,2) = T{27.0};
   A(2,0) = T{6.0}; A(2,1) = T{12.0}; A(2,2) = T{18.0};

   //A12..
   A(0,3) = T{8.0}; A(0,4) = T{16.0}; A(0,5) = T{24.0};
   A(1,3) = T{2.0}; A(1,4) = T{4.0}; A(1,5) = T{6.0};
   A(2,3) = T{3.0}; A(2,4) = T{6.0}; A(2,5) = T{9.0};

   //A13..
   A(0,6) = T{3.0}; A(0,7) = T{6.0}; A(0,8) = T{9.0};
   A(1,6) = T{0.0}; A(1,7) = T{0.0}; A(1,8) = T{0.0};
   A(2,6) = T{2.0}; A(2,7) = T{4.0}; A(2,8) = T{6.0};

   //A21..
   A(3,0) = T{4.0}; A(3,1) = T{6.0}; A(3,2) = T{8.0};
   A(4,0) = T{36.0}; A(4,1) = T{54.0}; A(4,2) = T{72.0};
   A(5,0) = T{24.0}; A(5,1) = T{36.0}; A(5,2) = T{48.0};

   //A22..
   A(3,3) = T{32.0}; A(3,4) = T{48.0}; A(3,5) = T{64.0};
   A(4,3) = T{8.0}; A(4,4) = T{12.0}; A(4,5) = T{16.0};
   A(5,3) = T{12.0}; A(5,4) = T{18.0}; A(5,5) = T{24.0};

   //A23..
   A(3,6) = T{12.0}; A(3,7) = T{18.0}; A(3,8) = T{24.0};
   A(4,6) = T{0.0}; A(4,7) = T{0.0}; A(4,8) = T{0.0};
   A(5,6) = T{8.0}; A(5,7) = T{12.0}; A(5,8) = T{16.0};

   //A31..
   A(6,0) = T{5.0}; A(6,1) = T{7.0}; A(6,2) = T{9.0};
   A(7,0) = T{45.0}; A(7,1) = T{63.0}; A(7,2) = T{81.0};
   A(8,0) = T{30.0}; A(8,1) = T{42.0}; A(8,2) = T{54.0};

   //A32..
   A(6,3) = T{40.0}; A(6,4) = T{56.0}; A(6,5) = T{72.0};
   A(7,3) = T{10.0}; A(7,4) = T{14.0}; A(7,5) = T{18.0};
   A(8,3) = T{15.0}; A(8,4) = T{21.0}; A(8,5) = T{27.0};

   //A33..
   A(6,6) = T{15.0}; A(6,7) = T{21.0}; A(6,8) = T{27.0};
   A(7,6) = T{0.0}; A(7,7) = T{0.0}; A(7,8) = T{0.0};
   A(8,6) = T{10.0}; A(8,7) = T{14.0}; A(8,8) = T{18.0};


   auto K = computeProductInf(G, R);
   T error_R = (A -  K).norm();

   if(error_R > prec){
      std::cout << "prod_inf, error P: " << error_R << std::endl;
      std::cout << "P analytique: " << A << std::endl;
      std::cout << "P calculé: " << K << std::endl;
   }

}


template< typename T>
void
test_identity(const T prec)
{
   static_matrix<T, 3, 3> I = static_matrix<T, 3, 3>::Identity();

   auto Isup = computeProductSup(I, I);
   T error_R = (Isup -  compute_IdentityTensor<T,3>()).norm();

   if(error_R > prec){
      std::cout << "Identity, error P: " << error_R << std::endl;
      std::cout << "P analytique: " << compute_IdentityTensor<T,3>() << std::endl;
      std::cout << "P calculé: " << Isup << std::endl;
   }

}
