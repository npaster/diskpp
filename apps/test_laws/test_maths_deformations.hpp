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
#include "BehaviorLaws/maths_deformation_tensors.hpp"

#pragma once


template< typename T>
void
test_Ftensor(const T prec)
{
   static_matrix<T, 3, 3> G;
   G(0,0) = T{1.0}; G(0,1) = T{2.0}; G(0,2) = T{3.0};
   G(1,0) = T{4.0}; G(1,1) = T{6.0}; G(1,2) = T{8.0};
   G(2,0) = T{5.0}; G(2,1) = T{7.0}; G(2,2) = T{9.0};
 
   static_matrix<T, 3, 3> F;
   F(0,0) = T{2.0}; F(0,1) = T{2.0}; F(0,2) = T{3.0};
   F(1,0) = T{4.0}; F(1,1) = T{7.0}; F(1,2) = T{8.0};
   F(2,0) = T{5.0}; F(2,1) = T{7.0}; F(2,2) = T{10.0};
   

   T error_R = (F - compute_FTensor(G)).norm();
   
   if(error_R > prec){
      std::cout << "error F computation: " << error_R << std::endl;
      std::cout << "F analytique: " << F << std::endl;
      std::cout << "F calculé: " << compute_FTensor(G) << std::endl;
   }
}


template< typename T>
void
test_CauchyGreen(const T prec)
{
   static_matrix<T, 3, 3> F;
   F(0,0) = T{2.0}; F(0,1) = T{2.0}; F(0,2) = T{3.0};
   F(1,0) = T{4.0}; F(1,1) = T{7.0}; F(1,2) = T{8.0};
   F(2,0) = T{5.0}; F(2,1) = T{7.0}; F(2,2) = T{10.0};
   
   static_matrix<T, 3, 3> C = F.transpose() * F;
   
   
   T error_R = (C - compute_CauchyGreenRightTensor(F)).norm();
   
   if(error_R > prec){
      std::cout << "error C computation: " << error_R << std::endl;
      std::cout << "C analytique: " << F << std::endl;
      std::cout << "C calculé: " << compute_CauchyGreenRightTensor(F) << std::endl;
   }
}