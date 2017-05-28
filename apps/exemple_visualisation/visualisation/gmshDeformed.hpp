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
#include <cassert>
#include <math.h>

#include "geometry/geometry.hpp"
#include "mesh/point.hpp"
#include "gmshMesh.hpp"
#include "gmshElement.hpp"
#include "gmshConvertMesh.hpp"

namespace visu{



template<typename T, size_t DIM>
void init_coor( const point<T,DIM>& point, std::vector<double>& coor)
{
   for(size_t i = 0; i < DIM; i++){
      coor[i] = double(point.at(i));
   }
}


template< size_t DIM>
void add_element( Gmesh& gmsh, const std::vector<Node>& list_nodes)
{
   if(DIM == 1) {
      assert(list_nodes.size() == 2);
      std::vector<size_t> index_nodes;
      for (size_t i = 0; i < list_nodes.size(); i++) {
         index_nodes.push_back(list_nodes[i].getIndex());
      }

      Edge edge(index_nodes, gmsh.getNumberofElements() +1, 0, 0 );
      gmsh.addEdge(edge);
   }
   else if(DIM == 2) {
      assert(list_nodes.size() >= 3);
      // sort node
      std::vector<size_t> nodes_sorted = sort_nodes_2D(list_nodes);

      std::vector<size_t> index_nodes;
      for (size_t i = 0; i < list_nodes.size(); i++) {
         index_nodes.push_back(list_nodes[nodes_sorted[i]].getIndex());
      }

      if(list_nodes.size() == 3) {
         Triangle tri(index_nodes, gmsh.getNumberofElements() +1, 0, 0 );
         gmsh.addTriangle(tri);
      }
      else if(list_nodes.size() == 4) {
         Quadrangle quad(index_nodes, gmsh.getNumberofElements() +1, 0, 0 );
         gmsh.addQuadrangle(quad);
      }
      else {
         std::vector<size_t> index_nodes_tri(3,0);
         index_nodes_tri[0] = index_nodes[0];

         for(size_t i = 1; i < (list_nodes.size()-1); i++){
            index_nodes_tri[1] = index_nodes[i];
            index_nodes_tri[2] = index_nodes[i+1];
            Triangle tri(index_nodes_tri, gmsh.getNumberofElements() +1, 0, 0 );
            gmsh.addTriangle(tri);
         }

      }
   }
   else if(DIM == 3) {
      assert(list_nodes.size() >= 4);
      // sort node
      std::vector<size_t> nodes_sorted = sort_nodes_3D(list_nodes);

      std::vector<size_t> index_nodes;
      for (size_t i = 0; i < list_nodes.size(); i++) {
         index_nodes.push_back(list_nodes[nodes_sorted[i]].getIndex());
      }

      if(list_nodes.size() == 4){// this is a tetrahedra
         Tetrahedron new_tetra(index_nodes, gmsh.getNumberofElements() +1, 0, 0 );
         gmsh.addTetrahedron(new_tetra);
      }
      else if(list_nodes.size() == 8){// this is a tetrahedra
         Hexahedron new_hexa(index_nodes, gmsh.getNumberofElements() +1, 0, 0 );
         gmsh.addHexahedron(new_hexa);
      }
   }
   else
      assert(false);
}


} //visu

#endif
