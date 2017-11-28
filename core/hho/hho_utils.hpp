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

namespace hho {

// scalar case
template<typename T>
T
eval(const dynamic_vector<T>& tab_coeff, const dynamic_vector<T>& base)
{
   assert(tab_coeff.rows() == base.rows());

   return tab_coeff.dot(base);
}

template<typename T, int DIM>
static_vector<T, DIM>
eval(const dynamic_vector<T>& tab_coeff, const Eigen::Matrix<T, Eigen::Dynamic, DIM>& base)
{
   assert(tab_coeff.rows() == base.rows());

   const auto prod = tab_coeff * base;

   static_vector<T, DIM> ret = static_vector<T, DIM>::Zero();

   for (size_t i = 0; i < DIM; i++) {
      ret(i) = prod(i);
   }

   return ret;
}

// vectorial case
template<typename T, int DIM>
static_matrix<T, DIM, DIM>
eval(const dynamic_vector<T>& tab_coeff, const std::vector<static_matrix<T, DIM, DIM>>& base)
{
   static_matrix<T, DIM, DIM> ret = static_matrix<T, DIM, DIM>::Zero();
   assert(tab_coeff.size() == base.size());

   for (size_t i = 0; i < base.size(); i++) {
      ret += tab_coeff(i) * base[i];
   }

   return ret;
}

// matricial case
template<typename T, int DIM>
static_vector<T, DIM>
eval(const dynamic_vector<T>& tab_coeff, const std::vector<static_vector<T, DIM>>& base)
{
   static_vector<T, DIM> ret = static_vector<T, DIM>::Zero();
   assert(tab_coeff.size() == base.size());

   for (size_t i = 0; i < base.size(); i++) {
      ret += tab_coeff(i) * base[i];
   }

   return ret;
}

} // hho namespace
} // disk namespace