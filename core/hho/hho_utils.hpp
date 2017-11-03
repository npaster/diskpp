/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

#pragma once

#include "common/eigen.hpp"

namespace disk {
   
   namespace hho{
      
      template<typename T, int DIM>
      static_matrix<T,DIM, DIM>
      eval_gradient(const dynamic_vector<T>& gradrec_coeff,
                    const std::vector<static_matrix<T,DIM, DIM> >& base_grad)
      {
         static_matrix<T,DIM, DIM> ret = static_matrix<T,DIM, DIM>::Zero();
         assert(gradrec_coeff.size() == base_grad.size());
         
         for(size_t i = 0; i < base_grad.size(); i++){
            ret += gradrec_coeff(i) * base_grad[i];
         }
         
         return ret;
      }
      
      template<typename T, int DIM>
      static_vector<T,DIM>
      eval_gradient(const dynamic_vector<T>& gradrec_coeff,
                    const std::vector<static_vector<T,DIM> >& base_grad)
      {
         static_vector<T,DIM> ret = static_vector<T,DIM>::Zero();
         assert(gradrec_coeff.size() == base_grad.size());
         
         for(size_t i = 0; i < base_grad.size(); i++){
            ret += gradrec_coeff(i) * base_grad[i];
         }
         
         return ret;
      }
      
   } // hho namespace
} // disk namespace