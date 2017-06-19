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
#include <regex>

#include "test_plaplace.hpp"
#include "test_tm_prod.hpp"
#include "test_maths_deformations.hpp"


int main(int argc, char **argv)
{
    using RealType = double;
    
    test_plaplace_2d<RealType>(1E-15);
    test_plaplace_3d<RealType>(1E-15);
    
    test_tm_prod<RealType>(1E-15);
    
    test_Ftensor<RealType>(1E-15);
    
    test_CauchyGreen<RealType>(1E-15);
    
    test_kronecker<RealType>(1E-15);
    test_prodsup<RealType>(1E-15);
    test_prodinf<RealType>(1E-15);
    test_identity<RealType>(1E-15);
}
