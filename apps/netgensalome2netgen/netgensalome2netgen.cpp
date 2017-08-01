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
#include <fstream>
#include <regex>
#include <array>
#include <vector>
#include "../../config.h"
#include "loaders/loader.hpp"

int main(int argc, char **argv)
{

    char *filename = argv[1];

    if (std::regex_match(filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: MEDIT 3D" << std::endl;

        std::ifstream   ifs(filename);
        std::string     keyword;
        
        if (!ifs.is_open())
        {
           std::cout << "Error opening " << filename << std::endl;
           return 0;
        }
        
        ifs >> keyword;
        if ( keyword != "MeshVersionFormatted" )
        {
           std::cout << "Expected keyword \"MeshVersionFormatted\"" << std::endl;
           return 0;
        }
        
        size_t format;
        ifs >> format;
        
        if ( format != 2 )
        {
           std::cout << "Expected format 2 (here: " << format << ")" << std::endl;
           return 0;
        }
        
        ifs >> keyword;
        if ( keyword != "Dimension" )
        {
           std::cout << "Expected keyword \"Dimension\"" << std::endl;
           return 0;
        }
        
        size_t dim;
        ifs >> dim;
        
        if (dim != 3 )
        {
           std::cout << "Expected dimension = 3 (here: " << dim << ")" << std::endl;
           return 0;
        }
        
        std::vector<std::array<double,3>> vertices;
        std::vector<std::array<size_t,4>> triangles;
        std::vector<std::array<size_t,5>> tetra;
        size_t elements_to_read(0);
        
        ifs >> keyword;
        while( keyword != "End")
        {
           std::cout << keyword  << std::endl;
           if ( keyword == "Vertices" )
           {
              ifs >> elements_to_read;
              vertices.reserve(elements_to_read);
              
              std::array<double,3> vertice = {0.0,0.0,0.0};
              
              for (size_t i = 0; i < elements_to_read; i++)
              {
                 ifs >> vertice[0] >> vertice[1] >> vertice[2] >> keyword;
                 vertices.push_back(vertice);
              }
                 
               
           }
           else if ( keyword == "Triangles" )
           {
              ifs >> elements_to_read;
              triangles.reserve(elements_to_read);
              
              std::array<size_t,4> tri = {0,0,0,0};
              
              for (size_t i = 0; i < elements_to_read; i++)
              {
                 ifs >> tri[0] >> tri[1] >> tri[2] >> tri[3];
                 triangles.push_back(tri);
              }
           }
           else if ( keyword == "Quadrilaterals" )
           {
              throw std::invalid_argument("mesh is not tetrahedral");
           }
           else if ( keyword == "Edges" )
           {
              ifs >> elements_to_read;
              
              std::array<size_t,3> edge = {0,0,0};
              
              for (size_t i = 0; i < elements_to_read; i++)
              {
                 ifs >> edge[0] >> edge[1] >> edge[2];
              }
           }
           else if ( keyword == "Tetrahedra" )
           {
              ifs >> elements_to_read;
              tetra.reserve(elements_to_read);
              
              std::array<size_t,5> tet = {0,0,0,0,0};
              
              for (size_t i = 0; i < elements_to_read; i++)
              {
                 ifs >> tet[0] >> tet[1] >> tet[2] >> tet[3] >> tet[4];
                 tetra.push_back(tet);
              }
           }
           else if ( keyword == "Hexahedra" )
           {
              throw std::invalid_argument("mesh is not tetrahedral");
           }
           else
           {
              std::cout << "Error parsing Medit file" << std::endl;
              return 0;
           }
           
           ifs >> keyword;
        }
        
        ifs.close();

        std::ofstream ofs(argv[2]);

        ofs << vertices.size() << std::endl;

        for (auto vertice : vertices)
           ofs << vertice[0] << " " << vertice[1] << " " << vertice[2] << std::endl;
        
        ofs << tetra.size() << std::endl;
        
        for (auto tet : tetra)
           ofs << tet[4] << " " << tet[0] << " " << tet[1] << " " << tet[2] << " " << tet[3] << std::endl;
        
        ofs << triangles.size() << std::endl;
        
        for (auto tri : triangles)
           ofs << tri[3] << " " << tri[0] << " " << tri[1] << " " << tri[2] << std::endl;


        ofs.close();
        
        return 1;
    }
    else
       std::cout << "Unknown format" << std::endl;
}
