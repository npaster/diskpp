/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

#pragma once

#include <vector>

namespace disk {

namespace mechanics {

/**
 * @brief Aitken acceleration
 *
 */
template < typename T >
class ConvergenceAcceleration {

    typedef dynamic_vector< T > vector_type;

    int n_iter;

    vector_type va_km, va_k;
    vector_type v_k, v_km;

  public:
    ConvergenceAcceleration() : n_iter( 0 ) {}

    vector_type aitken( const vector_type &v_kp ) {
        if ( n_iter == 0 ) {
            n_iter++;
            va_km = v_kp;
            return va_km;
        } else if ( n_iter == 1 ) {
            n_iter++;
            v_k = v_kp;
            va_k = v_kp;
            return v_k;
        } else {
            n_iter++;
            const vector_type va_d = va_k - va_km;
            const vector_type v_d = ( v_kp - va_k ) + ( v_k - va_km );

            // compute relaxation
            const T wr = va_d.dot( v_d ) / v_d.squaredNorm();

            // compute accelerted solution
            const vector_type va_kp = wr * v_kp + ( 1 - wr ) * va_k;

            // update
            v_k = v_kp;
            va_km = va_k;
            va_k = va_kp;

            return va_kp;
        }
    }

    vector_type aitken2( const vector_type &v_kp ) {
        if ( n_iter == 0 ) {
            n_iter++;
            v_km = v_kp;
            return v_km;
        } else if ( n_iter == 1 ) {
            n_iter++;
            v_k = v_kp;
            return v_k;
        } else {
            n_iter++;
            const vector_type vt = v_kp - 2 * v_k + v_km;
            const vector_type dv = v_kp - v_k;

            // compute relaxation
            const T wr = dv.dot( vt ) / vt.squaredNorm();

            // compute accelerted solution
            const vector_type va_kp = wr * v_k + ( 1 - wr ) * v_kp;

            // update
            v_km = v_k;
            v_k = v_kp;

            return va_kp;
        }
    }

    vector_type relaxation( const vector_type &v_kp, const T omega = 0.5 ) {
        if ( n_iter == 0 ) {
            n_iter++;
            va_k = v_kp;
            return va_k;
        } else {
            n_iter++;

            // compute accelerted solution
            const vector_type va_kp = ( 1.0 - omega ) * va_k + omega * v_kp;

            // update
            va_k = va_kp;

            return va_kp;
        }
    }
};
} // namespace mechanics
} // namespace disk