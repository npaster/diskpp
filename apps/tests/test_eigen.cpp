/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet (C) 2018         nicolas.pignet@enpc.fr
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

#include "common/eigen.hpp"
#include "verify.h"

int main()
{
    using T = double;

    // test eigen matrix
    static_matrix<T, 2, 3> mat;

    VERIFY_IS_EQUAL(mat.rows(), 2);
    VERIFY_IS_EQUAL(mat.cols(), 3);
    VERIFY_IS_EQUAL(mat.size(), 6);

    // test eigen vector
    static_vector<T, 5> vec;

    VERIFY_IS_EQUAL(vec.rows(), 5);
    VERIFY_IS_EQUAL(vec.cols(), 1);
    VERIFY_IS_EQUAL(vec.size(), 5);

    // test eigen 4th-tensor
    static_tensor<T, 3> tens;

    VERIFY_IS_EQUAL(tens.rows(), 9);
    VERIFY_IS_EQUAL(tens.cols(), 9);
    VERIFY_IS_EQUAL(tens.size(), 81);

    static_tensor<T, 2> tens2;

    VERIFY_IS_EQUAL(tens2.rows(), 4);
    VERIFY_IS_EQUAL(tens2.cols(), 4);
    VERIFY_IS_EQUAL(tens2.size(), 16);

    // test cross product
    static_vector<T, 2> v1;
    v1(0) = T(1.1);
    v1(1) = T(2.3);

    static_vector<T, 2> v2;
    v2(0) = T(2.1);
    v2(1) = T(1.3);

    static_vector<T, 3> v3 = cross(v1, v2);
    static_vector<T, 3> v4 = cross(v2, v1);

    VERIFY_IS_EQUAL(v3(0), T(0.0));
    VERIFY_IS_EQUAL(v3(1), T(0.0));
    VERIFY_IS_EQUAL(v3(2), -T(3.4));

    VERIFY_IS_EQUAL(v4(0), T(0.0));
    VERIFY_IS_EQUAL(v4(1), T(0.0));
    VERIFY_IS_EQUAL(v4(2), -v3(2));

    return 0;
}