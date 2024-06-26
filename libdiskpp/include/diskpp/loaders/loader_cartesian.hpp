/*
 *       /\        Nicolas Pignet (C) 2024
 *      /__\       nicolas.pignet@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

#pragma once

#include <fstream>
#include <thread>
#include <tuple>
#include <vector>

#include "diskpp/common/mapped_file.h"
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/loaders/mesh_loader.hpp"
#include "diskpp/loaders/strtot.hpp"
#include "diskpp/mesh/mesh.hpp"

namespace disk
{

namespace priv
{

template<typename T>
std::tuple<T, T, T, T, T, T, T, T>
read_hexahedron_line(const char* str, char** endptr)
{
    T t1, t2, t3, t4, t5, t6, t7, t8;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);
    t5 = strtot<T>(*endptr, endptr);
    t6 = strtot<T>(*endptr, endptr);
    t7 = strtot<T>(*endptr, endptr);
    t8 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2, t3, t4, t5, t6, t7, t8);
}

template<typename T>
std::tuple<T, T, T, T>
read_hex_face_line(const char* str, char** endptr)
{
    T t1, t2, t3, t4;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2, t3, t4);
}

template<typename T>
std::tuple<T, T, T, T>
read_quad_line(const char* str, char** endptr)
{
    T t1, t2, t3, t4;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2, t3, t4);
}

template<typename T>
std::tuple<T, T>
read_quad_face_line(const char* str, char** endptr)
{
    T t1, t2;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2);
}

} // namespace priv

template<typename T, size_t DIM>
class cartesian_mesh_loader;

template<typename T>
class cartesian_mesh_loader<T, 3> : public mesh_loader<cartesian_mesh<T, 3>>
{
    typedef cartesian_mesh<T, 3>             mesh_type;
    typedef typename mesh_type::point_type   point_type;
    typedef typename mesh_type::node_type    node_type;
    typedef typename mesh_type::edge_type    edge_type;
    typedef typename mesh_type::surface_type surface_type;
    typedef typename mesh_type::volume_type  volume_type;

    std::vector<point_type>   points;
    std::vector<node_type>    nodes;
    std::vector<edge_type>    edges;
    std::vector<surface_type> surfaces, boundary_surfaces;
    std::vector<volume_type>  volumes;

    bool
    hex_read(const std::string& filename)
    {
        /* Open file */
        if (filename.size() == 0)
        {
            std::cout << "Invalid mesh file name" << std::endl;
            return false;
        }

        size_t lines, linecount;

        mapped_file mf(filename);

        // std::cout << green << " * * * Reading NETGEN format mesh * * * ";
        // std::cout << nocolor << std::endl;

        /************************ Read points ************************/
        linecount = 0;

        const char* data = mf.mem();
        char*       endptr;

        lines = strtot<size_t>(data, &endptr);

        points.reserve(lines);
        nodes.reserve(lines);

        while (linecount < lines)
        {
            if (this->verbose() && ((linecount % 100000) == 0))
            {
                std::cout << "Reading points: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto point = priv::read_3d_point_line<T>(endptr, &endptr, 1.0);

            points.push_back(point);

            auto point_id = disk::point_identifier<3>(linecount);
            auto node     = node_type({point_id});

            nodes.push_back(node);
            /* Do something with that point */

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading points: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        /************************ Read hexahedra ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        edges.reserve(lines * 12);
        surfaces.reserve(lines * 6);
        volumes.reserve(lines);

        while (linecount < lines)
        {
            if (this->verbose() && ((linecount % 100000) == 0))
            {
                std::cout << "Reading hexahedra: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_hexahedron_line<size_t>(endptr, &endptr);

            disk::point_identifier<3> p0(std::get<0>(t));
            disk::point_identifier<3> p1(std::get<1>(t));
            disk::point_identifier<3> p2(std::get<2>(t));
            disk::point_identifier<3> p3(std::get<3>(t));
            disk::point_identifier<3> p4(std::get<4>(t));
            disk::point_identifier<3> p5(std::get<5>(t));
            disk::point_identifier<3> p6(std::get<6>(t));
            disk::point_identifier<3> p7(std::get<7>(t));

            edges.push_back(edge_type({p0, p1}));
            edges.push_back(edge_type({p0, p2}));
            edges.push_back(edge_type({p0, p4}));
            edges.push_back(edge_type({p1, p3}));
            edges.push_back(edge_type({p1, p5}));
            edges.push_back(edge_type({p2, p3}));
            edges.push_back(edge_type({p2, p6}));
            edges.push_back(edge_type({p3, p7}));
            edges.push_back(edge_type({p4, p5}));
            edges.push_back(edge_type({p4, p6}));
            edges.push_back(edge_type({p5, p7}));
            edges.push_back(edge_type({p6, p7}));

            surfaces.push_back(surface_type({p0, p2, p6, p4}));
            surfaces.push_back(surface_type({p1, p3, p7, p5}));
            surfaces.push_back(surface_type({p0, p1, p3, p2}));
            surfaces.push_back(surface_type({p4, p5, p7, p6}));
            surfaces.push_back(surface_type({p0, p4, p5, p1}));
            surfaces.push_back(surface_type({p2, p6, p7, p3}));

            volumes.push_back(volume_type({p0, p1, p2, p3, p4, p5, p6, p7}));

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading hexahedra: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        /************************ Read boundary surfaces ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        boundary_surfaces.reserve(lines);

        while (linecount < lines)
        {
            if (this->verbose() && ((linecount % 50000) == 0))
            {
                std::cout << "Reading hex face: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_hex_face_line<size_t>(endptr, &endptr);

            disk::point_identifier<3> p0(std::get<0>(t));
            disk::point_identifier<3> p1(std::get<1>(t));
            disk::point_identifier<3> p2(std::get<2>(t));
            disk::point_identifier<3> p3(std::get<3>(t));

            surface_type quad({p0, p1, p2, p3});

            boundary_surfaces.push_back(quad);

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading hex face: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        return true;
    }

  public:
    static const char constexpr* expected_extension = "hex";
    cartesian_mesh_loader()                         = default;

    bool
    read_mesh(const std::string& s)
    {
        if (this->verbose())
            std::cout << " *** READING CARTESIAN 3D MESH ***" << std::endl;

        return hex_read(s);
    }

    bool
    populate_mesh(mesh_type& msh)
    {
        auto storage = msh.backend_storage();

        if (this->verbose())
        {
            std::cout << "Sorting data...";
            std::cout.flush();
        }

        storage->points = std::move(points);
        storage->nodes  = std::move(nodes);

        /* sort edges, make unique and move them in geometry */
        THREAD(edge_thread, priv::sort_uniq(edges); storage->edges = std::move(edges););

        /* sort triangles, make unique and move them in geometry */
        THREAD(quad_thread, priv::sort_uniq(surfaces); storage->surfaces = std::move(surfaces););

        /* sort tetrahedra, make unique and move them in geometry */
        THREAD(hex_thread, std::sort(volumes.begin(), volumes.end()); storage->volumes = std::move(volumes););

        /* wait for the threads */
        WAIT_THREAD(edge_thread);
        WAIT_THREAD(quad_thread);
        WAIT_THREAD(hex_thread);

        storage->boundary_info.resize(storage->surfaces.size());
        for (auto& bs : boundary_surfaces)
        {
            auto position = find_element_id(storage->surfaces.begin(), storage->surfaces.end(), bs);
            if (position.first == false)
            {
                std::cout << "Bad bug at " << __FILE__ << "(" << __LINE__ << ")" << std::endl;
                return false;
            }

            boundary_descriptor bi(0, true);
            storage->boundary_info.at(position.second) = bi;
        }

        if (this->verbose())
        {
            std::cout << "done." << std::endl;

            std::cout << "Nodes: " << storage->nodes.size() << std::endl;
            std::cout << "Edges: " << storage->edges.size() << std::endl;
            std::cout << "Faces: " << storage->surfaces.size() << std::endl;
            std::cout << "Volumes: " << storage->volumes.size() << std::endl;
        }

        boundary_surfaces.clear();

        return true;
    }
};

template<typename T>
class cartesian_mesh_loader<T, 2> : public mesh_loader<cartesian_mesh<T, 2>>
{
    typedef cartesian_mesh<T, 2>             mesh_type;
    typedef typename mesh_type::point_type   point_type;
    typedef typename mesh_type::node_type    node_type;
    typedef typename mesh_type::edge_type    edge_type;
    typedef typename mesh_type::surface_type surface_type;

    std::vector<point_type>   points;
    std::vector<node_type>    nodes;
    std::vector<edge_type>    edges, boundary_edges;
    std::vector<surface_type> surfaces;

    bool
    quad_read(const std::string& filename)
    {
        /* Open file */
        if (filename.size() == 0)
        {
            std::cout << "Can't open '" << filename << "'" << std::endl;
            return false;
        }

        size_t lines, linecount;

        mapped_file mf(filename);
        if (!mf.is_open())
            return false;

        /************************ Read points ************************/
        linecount = 0;

        const char* data = mf.mem();
        char*       endptr;

        lines = strtot<size_t>(data, &endptr);

        points.reserve(lines);
        nodes.reserve(lines);

        while (linecount < lines)
        {
            if (this->verbose() && (linecount % 100000) == 0)
            {
                std::cout << "Reading points: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto point = priv::read_2d_point_line<T>(endptr, &endptr, 1.0);

            points.push_back(point);

            auto point_id = disk::point_identifier<2>(linecount);
            auto node     = node_type({point_id});

            nodes.push_back(node);
            /* Do something with that point */

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading points: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        /************************ Read hexahedra ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        edges.reserve(lines * 4);
        surfaces.reserve(lines);

        while (linecount < lines)
        {
            if (this->verbose() && (linecount % 100000) == 0)
            {
                std::cout << "Reading quads: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_quad_line<size_t>(endptr, &endptr);

            disk::point_identifier<2> p0(std::get<0>(t));
            disk::point_identifier<2> p1(std::get<1>(t));
            disk::point_identifier<2> p2(std::get<2>(t));
            disk::point_identifier<2> p3(std::get<3>(t));

            edges.push_back(edge_type({p0, p1}));
            edges.push_back(edge_type({p0, p2}));
            edges.push_back(edge_type({p1, p3}));
            edges.push_back(edge_type({p2, p3}));

            surfaces.push_back(surface_type({p0, p1, p2, p3}));

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading quads: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        /************************ Read boundary surfaces ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        boundary_edges.reserve(lines);

        while (linecount < lines)
        {
            if (this->verbose() && (linecount % 50000) == 0)
            {
                std::cout << "Reading faces: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_quad_face_line<size_t>(endptr, &endptr);

            disk::point_identifier<2> p0(std::get<0>(t));
            disk::point_identifier<2> p1(std::get<1>(t));

            edge_type bnd({p0, p1});

            boundary_edges.push_back(bnd);

            linecount++;
        }

        if (this->verbose())
        {
            std::cout << "Reading faces: " << linecount;
            std::cout << "/" << lines << std::endl;
        }

        return true;
    }

  public:
    static const char constexpr* expected_extension = "quad";

    cartesian_mesh_loader() = default;

    bool
    read_mesh(const std::string& s)
    {
        if (this->verbose())
            std::cout << " *** READING CARTESIAN 2D MESH ***" << std::endl;

        return quad_read(s);
    }

    bool
    populate_mesh(mesh_type& msh)
    {
        auto storage = msh.backend_storage();

        if (this->verbose())
        {
            std::cout << "Sorting data...";
            std::cout.flush();
        }

        storage->points = std::move(points);
        storage->nodes  = std::move(nodes);

        /* sort edges, make unique and move them in geometry */
        THREAD(edge_thread, priv::sort_uniq(edges); storage->edges = std::move(edges););

        /* sort triangles, make unique and move them in geometry */
        THREAD(quad_thread, priv::sort_uniq(surfaces); storage->surfaces = std::move(surfaces););

        /* wait for the threads */
        WAIT_THREAD(edge_thread);
        WAIT_THREAD(quad_thread);

        storage->boundary_info.resize(storage->edges.size());
        for (auto& bs : boundary_edges)
        {
            auto position = find_element_id(storage->edges.begin(), storage->edges.end(), bs);
            if (position.first == false)
            {
                std::cout << "Bad bug at " << __FILE__ << "(" << __LINE__ << ")" << std::endl;
                return false;
            }

            boundary_descriptor bi(0, true);
            storage->boundary_info.at(position.second) = bi;
        }

        if (this->verbose())
        {
            std::cout << "done." << std::endl;

            std::cout << "Nodes: " << storage->nodes.size() << std::endl;
            std::cout << "Edges: " << storage->edges.size() << std::endl;
            std::cout << "Faces: " << storage->surfaces.size() << std::endl;
        }

        boundary_edges.clear();

        return true;
    }
};

} // end namespace disk