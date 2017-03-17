/*
 *       /\
 *      /__\       Nicolas Pignet (C) 2017
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    Gmsh tools
 *  /_\/_\/_\/_\
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#ifndef GMSHDISK_H
#define GMSHDISK_H

#include<vector>
#include<string>
#include<utility>
#include <array>
#include <cassert>

#include "loaders/loader.hpp"
#include "geometry/geometry.hpp"
#include "mesh/point.hpp"
#include "gmshMesh.h"
#include "gmshElement.h"

namespace visu{

template<typename T, size_t DIM>
Node convertPoint( const point<T,DIM>& point, const size_t index)
{
   std::vector<double> coor(3,0.0);

   for(size_t i = 0; i < DIM; i++){
      coor[i] = double(point.at(i));
   }

   Node node(coor, index, 0);

   return node;
}

template<typename DiskEdge>
Edge convertEdge( const DiskEdge& edge, const size_t index)
{
   auto nodes = edge.point_ids();

   assert(nodes.size() == 2);

   std::vector<size_t> tmpnodes = { nodes[0] + 1, nodes[1] + 1 };

   Edge tmpedge(tmpnodes, index, 0, 0);

   return tmpedge;
}

template<typename DiskTriangle>
Triangle convertTriangle( const DiskTriangle& triangle, const size_t index)
{
   auto nodes = triangle.point_ids();

   assert(nodes.size() == 3);

//    std::cout << " ------------------------------------------------ " << std::endl;
//    std::cout << nodes[0] + 1 << " " << nodes[1] + 1 << " "<< nodes[2] + 1 << " " << std::endl;

   std::vector<size_t> tmpnodes = { nodes[0] + 1, nodes[1] + 1, nodes[2] + 1 };

   Triangle tri(tmpnodes, index, 0, 0);

   return tri;
}

template<typename DiskQuadrangle>
Quadrangle convertQuadrangle( const DiskQuadrangle& quadrangle, const size_t index)
{
   auto nodes = quadrangle.point_ids();

   assert(nodes.size() == 4);

   std::vector<size_t> tmpnodes = { nodes[0] + 1, nodes[1] + 1, nodes[3] + 1, nodes[2] + 1 }; //the ordering of nodes a bit strange

   Quadrangle quad(tmpnodes, index, 0, 0);

   return quad;
}

template<typename DiskTetra>
Tetrahedron convertTetrahedron( const DiskTetra& tetrahedron, const size_t index)
{
   auto nodes = tetrahedron.point_ids();

   assert(nodes.size() == 4);
   /*
   std::cout << " ------------------------------------------------ " << std::endl;
   std::cout << nodes[0] + 1 << " " << nodes[1] + 1 << " "<< nodes[2] + 1 << " "<< nodes[3] + 1 << " " << std::endl;*/

   std::vector<size_t> tmpnodes = { nodes[0] + 1, nodes[1] + 1, nodes[2] + 1, nodes[3] + 1 };

   Tetrahedron tetra(tmpnodes, index, 0, 0);

   return tetra;
}


template<typename DiskHexa>
Hexahedron convertHexahedron( const DiskHexa& hexahedron, const size_t index)
{
   auto nodes = hexahedron.point_ids();

   assert(nodes.size() == 8);

   /*std::cout << " ------------------------------------------------ " << std::endl;
   std::cout << nodes[0] + 1 << " " << nodes[1] + 1 << " "<< nodes[3] + 1 << " "<< nodes[2] + 1 << " " << std::endl;
   std::cout << nodes[4] + 1 << " " << nodes[5] + 1 << " "<< nodes[7] + 1 << " "<< nodes[6] + 1 << " " << std::endl;*/

   std::vector<size_t> tmpnodes = { nodes[0] + 1, nodes[1] + 1, nodes[3] + 1, nodes[2] + 1,
                                    nodes[4] + 1, nodes[5] + 1, nodes[7] + 1, nodes[6] + 1 }; //the ordering of nodes a bit strange

   Hexahedron hexa(tmpnodes, index, 0, 0);

   return hexa;
}


} //visu

#endif
