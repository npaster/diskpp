/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016, 2017 - matteo.cicuttin@enpc.fr
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
 * cite
*/

#pragma once



enum BoundaryType : size_t
{
   CLAMPED = 0,
   DX = 1,
   DY = 2,
   DZ = 3,
   DXDY = 4,
   DXDZ = 5,
   DYDZ = 6
};


struct BoundaryConditions
{
   size_t id;
   size_t boundary_type;
   // add a function for each boundary and logical to say if we use it
};
   
   





