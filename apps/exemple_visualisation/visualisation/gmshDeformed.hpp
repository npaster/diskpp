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

#ifndef GMSHDEFORMED_H
#define GMSHDEFORMED_H

#include<vector>
#include<string>
#include <cassert>

#include "geometry/geometry.hpp"
#include "mesh/point.hpp"
#include "gmshMesh.h"
#include "gmshElement.h"
#include "gmshDisk.h"

namespace visu{



template<typename T, size_t DIM>
void init_coor( const point<T,DIM>& point, std::vector<double>& coor)
{
   for(size_t i = 0; i < DIM; i++){
      coor[i] = double(point.at(i));
   }
}

template< size_t DIM>
void add_element( Gmesh gmsh, const std::vector<Node> list_nodes)
{
   if(DIM == 1) {
      assert(list_nodes.size() == 2);
      std::vector<size_t> index_nodes;
      for (size_t i = 0; i < list_nodes.size(); i++) {
         gmsh.addNode(list_nodes[i]);
         index_nodes.push_back(list_nodes[i].index());
      }

      Edge edge(index_nodes, gmsh.getNumberofElements() +1, 0, 0 );
      gmsh.addEdge(edge);
   }
   else if(DIM == 2) {
      assert(list_nodes.size() >= 3);
      std::vector<size_t> index_nodes;
      for (size_t i = 0; i < list_nodes.size(); i++) {
         gmsh.addNode(list_nodes[i]);
         index_nodes.push_back(list_nodes[i].index());
      }

      //ajouter test orientation
      Triangle tri(index_nodes, gmsh.getNumberofElements() +1, 0, 0 );
      gmsh.addTriangle(tri);
   }
   else
      assert(false);
}


} //visu

#endif
