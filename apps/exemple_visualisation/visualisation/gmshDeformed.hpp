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

namespace visu{



template<typename T, size_t DIM>
void init_coor( const point<T,DIM>& point, std::vector<double>& coor)
{
   for(size_t i = 0; i < DIM; i++){
      coor[i] = double(point.at(i));
   }
}

std::vector<size_t>
sort_nodes(const std::vector<Node>& list_nodes)
{
   //Works only for 2D case
   //compute barycenter
   std::vector<double> bar(2,0.0);

   for(size_t i = 0; i < list_nodes.size(); i++) {
      std::vector<double> coor_node = list_nodes[i].getCoordinate();
      bar[0] += coor_node[0];
      bar[1] += coor_node[1];
   }

   bar[0] /= double(list_nodes.size());
   bar[1] /= double(list_nodes.size());

   //compute angle
   std::vector<double> angle(list_nodes.size(), 0.0);
   for(size_t i = 0; i < list_nodes.size(); i++) {
      std::vector<double> coor_node = list_nodes[i].getCoordinate();
      double x = coor_node[0] - bar[0];
      double y = coor_node[1] - bar[1];
      angle[i] = atan2(y,x);
   }

   std::vector<size_t> ret(list_nodes.size(), 0);

   for(size_t i = 0; i < list_nodes.size(); i++) {
      ret[i] = i;
   }

   for(size_t i = 0; i < (list_nodes.size()-1); i++) {
      double vmin = angle[ret[i]];
      size_t imin = i;

      for(size_t j = i+1; j < list_nodes.size(); j++) {
         if( angle[ret[j]] < vmin) {
            vmin = angle[ret[j]];
            imin = j;
         }
      }

      size_t itmp = ret[i];
      ret[i] = ret[imin];
      ret[imin] = itmp;
   }

   return ret;
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
      std::vector<size_t> nodes_sorted = sort_nodes(list_nodes);

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
         index_nodes_tri[1] = index_nodes[1];
         for(size_t i = 2; i < list_nodes.size(); i++){
            index_nodes_tri[2] = index_nodes[i];
            Triangle tri(index_nodes_tri, gmsh.getNumberofElements() +1, 0, 0 );
            gmsh.addTriangle(tri);
         }

      }
   }
   else
      assert(false);
}


} //visu

#endif
