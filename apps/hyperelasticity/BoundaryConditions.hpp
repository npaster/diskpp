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

enum DirichletType : size_t
{
   CLAMPED = 0,
   DX = 1,
   DY = 2,
   DZ = 3,
   DXDY = 4,
   DXDZ = 5,
   DYDZ = 6,
   OTHER = 7
};


enum NeumannType : size_t
{
   FREE = 0,
   NEUMANN = 1
};



struct BoundaryType
{
   size_t id;
   size_t boundary_type;
   // add a function for each boundary and logical to say if we use it
};


class BoundaryConditions
{
private:
   std::vector<BoundaryType>  m_neumann_conditions;
   std::vector<BoundaryType>  m_dirichlet_conditions;
   std::vector<std::pair<bool, size_t>>  m_faces_dirichlet;
   std::vector<std::array<size_t, 2>> m_lagranges_info;

   size_t   m_nb_faces_boundary;
   size_t   m_nb_faces_dirichlet;
   size_t   m_nb_faces_neumann;
   size_t   m_nb_lag;



   template<typename TypeMesh>
   void
   find_neumann_faces(const TypeMesh& msh)
   {
      m_faces_dirichlet.reserve(m_nb_faces_boundary);
      m_faces_dirichlet.assign(m_nb_faces_boundary, std::make_pair(true, OTHER) );
      m_nb_faces_dirichlet = m_nb_faces_boundary;
      m_nb_faces_neumann = 0;

      if(!m_neumann_conditions.empty()){
         size_t face(0);
         for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
         {
            auto bfc = *itor;

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first)
               throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;
            const size_t b_id = msh.boundary_id(face_id);

            //Find if this face is a boundary face with Neumann Condition
            for(auto& elem : m_neumann_conditions)
            {
               if(b_id == elem.id){
                  m_faces_dirichlet[face] = std::make_pair(false, elem.boundary_type);
                  m_nb_faces_dirichlet--;
                  m_nb_faces_neumann++;
                  break;
               }
            }
            face++;
         }
      }
   }

   template<typename TypeMesh>
   void
   number_of_lag_conditions(const TypeMesh& msh)
   {
      //By default all boundary faces are dirichlet faces
      const size_t DIM = msh.dimension;
      m_nb_lag = 0;

      if(!m_dirichlet_conditions.empty()){
         size_t face(0);
         std::array<size_t, 2> info = {0,0};
         for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
         {
            if(m_faces_dirichlet[face].first){
               auto bfc = *itor;

               auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
               if (!eid.first)
               throw std::invalid_argument("This is a bug: face not found");

               auto face_id = eid.second;
               bool dirichlet_standart = true;

               for(auto& elem : m_dirichlet_conditions)
               {
                  if(msh.boundary_id(face_id) == elem.id)
                  {
                     dirichlet_standart = false;
                     switch ( elem.boundary_type ) {
                        case CLAMPED:
                           info[0] = m_nb_lag; info[1] = DIM;
                           m_lagranges_info.push_back(info);
                           m_nb_lag += DIM;
                           m_faces_dirichlet[face].second = CLAMPED;
                           break;
                        case DX:
                           info[0] = m_nb_lag; info[1] = 1;
                           m_lagranges_info.push_back(info);
                           m_nb_lag += 1;
                           m_faces_dirichlet[face].second = DX;
                           break;
                        case DY:
                           info[0] = m_nb_lag; info[1] = 1;
                           m_lagranges_info.push_back(info);
                           m_nb_lag += 1;
                           m_faces_dirichlet[face].second = DY;
                           break;
                        case DZ:
                           if( DIM != 3){
                              std::cout << "Invalid condition for face:" << face_id << std::endl;
                              throw std::invalid_argument(" ONLY DIM = 3 for this Dirichlet Conditions");
                           }
                           else{
                              info[0] = m_nb_lag; info[1] = 1;
                              m_lagranges_info.push_back(info);
                              m_nb_lag += 1;
                              m_faces_dirichlet[face].second = DZ;
                           }
                           break;
                        case DXDY:
                           if( DIM != 3){
                              std::cout << "Invalid condition for face:" << face_id << std::endl;
                              throw std::invalid_argument(" ONLY DIM = 3 for this Dirichlet Conditions");
                           }
                           else{
                              info[0] = m_nb_lag; info[1] = 2;
                              m_lagranges_info.push_back(info);
                              m_nb_lag += 2;
                              m_faces_dirichlet[face].second = DXDY;
                           }
                           break;
                        case DXDZ:
                           if( DIM != 3){
                              std::cout << "Invalid condition for face:" << face_id << std::endl;
                              throw std::invalid_argument(" ONLY DIM = 3 for this Dirichlet Conditions");
                           }
                           else{
                              info[0] = m_nb_lag; info[1] = 2;
                              m_lagranges_info.push_back(info);
                              m_nb_lag += 2;
                              m_faces_dirichlet[face].second = DXDZ;
                           }
                           break;
                        case DYDZ:
                           if( DIM != 3){
                              std::cout << "Invalid condition for face:" << face_id << std::endl;
                              throw std::invalid_argument(" ONLY DIM = 3 for this Dirichlet Conditions");
                           }
                           else{
                              info[0] = m_nb_lag; info[1] = 2;
                              m_lagranges_info.push_back(info);
                              m_nb_lag += 2;
                              m_faces_dirichlet[face].second = DYDZ;
                           }
                           break;
                        default:
                           std::cout << "Unknown Dirichlet Conditions: we do nothing" << std::endl;
                           break;
                     }
                     break;
                  }
               }

               if(dirichlet_standart){
                  info[0] = m_nb_lag; info[1] = DIM;
                  m_lagranges_info.push_back(info);
                  m_nb_lag += DIM;
               }
            }
            face++;
         }
      }
   }

public:

   BoundaryConditions() : m_nb_faces_boundary(0), m_nb_faces_dirichlet(0),
                         m_nb_faces_neumann(0), m_nb_lag(0) {}

   template< typename MeshType>
   BoundaryConditions(const MeshType& msh,
                     const std::vector<BoundaryType>& neumann_conditions,
                     const std::vector<BoundaryType>& dirichlet_conditions) :
                     m_nb_faces_boundary(msh.boundary_faces_size()),
                     m_neumann_conditions(neumann_conditions),
                     m_dirichlet_conditions(dirichlet_conditions)
   {
      find_neumann_faces(msh);
      number_of_lag_conditions(msh);

      std::cout << "face: "  << m_nb_faces_boundary << '\n';
      std::cout << "dir: "  << m_nb_faces_dirichlet << '\n';
      std::cout << "neu: "  << m_nb_faces_neumann << '\n';
      std::cout << "lag: "  << m_nb_lag << '\n';
   }

   size_t nb_faces_boundary() const {return m_nb_faces_boundary;}
   size_t nb_faces_dirichlet() const {return m_nb_faces_dirichlet;}
   size_t nb_faces_neumann() const {return m_nb_faces_neumann;}
   size_t nb_lag() const {return m_nb_lag;}

   std::vector<BoundaryType>  boundary_neumann() const
   {
      return m_neumann_conditions;
   }

   bool is_boundary_dirichlet(const size_t face_i) const
   {
      return m_faces_dirichlet[face_i].first;
   }

   size_t boundary_type(const size_t face_i) const
   {
      return m_faces_dirichlet[face_i].second;
   }

   size_t nb_lag_conditions_faceI(const size_t face_i) const
   {
      return m_lagranges_info[face_i][1];
   }

   size_t begin_lag_conditions_faceI(const size_t face_i) const
   {
      return m_lagranges_info[face_i][0];
   }
};
