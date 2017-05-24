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

#pragma once


#include<vector>
#include<string>
#include<utility>
#include <array>
#include <cassert>

#include "loaders/loader.hpp"
#include "geometry/geometry.hpp"
#include "mesh/point.hpp"
#include "gmshMesh.hpp"
#include "gmshElement.hpp"
#include "gmshData.hpp"
#include "gmshDisk.hpp"

namespace visu{
   std::vector<size_t>
   sort_nodes_2D(const std::vector<Node>& list_nodes)
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

   template<typename Mesh>
   Gmesh convertMesh(const Mesh& mesh)
   {

      Gmesh msh;

      auto storage = mesh.backend_storage();

      // Find the dimension of the mesh
      size_t DIM = mesh.dimension;

      // Fill the gmsh's mesh
      // conversion point in node
      size_t nb_node(0);
      for(auto point :  storage->points){
         nb_node +=1;
         msh.addNode(convertPoint(point, nb_node));

         Vertice vert(nb_node, nb_node, 0, 0);
         msh.addVertice(vert);
      }
      assert(storage->points.size() == msh.getNumberofNodes());


      // conversion edge
      size_t nb_edge(0);
      for(auto edge :  storage->edges){
         nb_edge +=1;
         auto list_nodes = edge.point_ids();
         assert(list_nodes.size() == 2);
         std::vector<size_t> index_nodes;
         for (size_t i = 0; i < 2; i++) {
            index_nodes.push_back(list_nodes[i] + 1);
         }

         Edge new_edge(index_nodes, msh.getNumberofElements() +1, 0, 0 );
         msh.addEdge(new_edge);
      }
      assert(storage->edges.size() == msh.getEdges().size());

      if(DIM == 2){
         // conversion surface
         size_t nb_surfaces(0);
         for(auto surface :  storage->surfaces){
            nb_surfaces +=1;
            auto list_nodes_index = surface.point_ids();
            std::vector<Node> list_nodes;

            //recup the coordinate of nodes
            for (size_t i = 0; i < list_nodes_index.size(); i++) {
               auto pt = storage->points[list_nodes_index[i]];
               std::vector<double> coor(3, double{0});

               for (size_t i = 0; i < 2; i++) {
                  coor[i] = double(pt.at(i));
               }

               Node tmp_node(coor, i+1, 0);
               list_nodes.push_back(tmp_node);
            }


            std::vector<size_t> nodes_sorted = sort_nodes_2D(list_nodes);
            std::vector<size_t> index_nodes;
            for (size_t i = 0; i < list_nodes.size(); i++) {
               index_nodes.push_back(list_nodes_index[nodes_sorted[i]] + 1);
            }

            if(list_nodes.size() == 3) {
               Triangle new_tri(index_nodes, msh.getNumberofElements() +1, 0, 0 );
               msh.addTriangle(new_tri);
            }
            else if(list_nodes.size() == 4) {
               Quadrangle new_quad(index_nodes, msh.getNumberofElements() +1, 0, 0 );
               msh.addQuadrangle(new_quad);
            }
            else {
               std::vector<size_t> index_nodes_tri(3,0);
               index_nodes_tri[0] = index_nodes[0];

               for(size_t i = 1; i < (list_nodes.size()-1); i++){
                  index_nodes_tri[1] = index_nodes[i];
                  index_nodes_tri[2] = index_nodes[i+1];
                  Triangle new_tri(index_nodes_tri, msh.getNumberofElements() +1, 0, 0 );
                  msh.addTriangle(new_tri);
               }

            }
         }
         assert(storage->surfaces.size() == nb_surfaces);
      }

      return msh;

   };
} //visu
