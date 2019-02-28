/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
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

#include <iostream>
#include <sstream>

#include <vector>

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "bases/bases.hpp"

namespace disk
{

template<typename MeshType>
class gmsh_io
{
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    const static size_t dimension = mesh_type::dimension;

    PostMesh<mesh_type> m_post_msh;
    gmsh::Gmesh         m_gmsh;

  public:
    gmsh_io() {}

    gmsh_io(const mesh_type& msh)
    {
        // compute simplicial submesh for post-processing
        m_post_msh = disk::PostMesh<mesh_type>(msh);
        m_gmsh     = disk::convertMesh(m_post_msh);
    }

    PostMesh<mesh_type>&
    post_mesh()
    {
        return m_post_msh;
    }

    PostMesh<mesh_type>
    post_mesh() const
    {
        return m_post_msh;
    }

    gmsh::Gmesh&
    gmesh()
    {
        return m_gmsh;
    }

    gmsh::Gmesh
    gmesh() const
    {
        return m_gmsh;
    }

    void
    save_mesh(const std::string& filename) const
    {
        // Save mesh
        m_gmsh.writeGmesh(filename, 2);
    }

    template<typename Function>
    void
    plot_scalar_function(const std::string& filename, const Function& fct) const
    {
        std::vector<gmsh::Data>    data;         // create data
        std::vector<gmsh::SubData> subdata;      // create subdata
        data.reserve(m_gmsh.getNumberofNodes()); // data has a size of nb_node

        auto storage = m_post_msh.mesh().backend_storage();

        // Compute the average value and save it
        for (int i_node = 1; i_node <= m_gmsh.getNumberofNodes(); i_node++)
        {
            const auto pt = storage->points[i_node - 1];

            const auto sol = fct(pt);

            const gmsh::Data tmp_data(i_node, disk::convertToVectorGmsh(sol));
            data.push_back(tmp_data);
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(1, 0.0, "sol_node", data, subdata);
        // Save the view
        nodedata.saveNodeData(filename, m_gmsh);
    }

    template<typename Function>
    void
    plot_vector_function(const std::string& filename, const Function& fct) const
    {
        std::vector<gmsh::Data>    data;         // create data
        std::vector<gmsh::SubData> subdata;      // create subdata
        data.reserve(m_gmsh.getNumberofNodes()); // data has a size of nb_node

        auto storage = m_post_msh.mesh().backend_storage();

        // Compute the average value and save it
        for (int i_node = 1; i_node <= m_gmsh.getNumberofNodes(); i_node++)
        {
            const auto pt = storage->points[i_node - 1];

            const auto sol  = fct(pt);

            const gmsh::Data tmp_data(i_node, disk::convertToVectorGmsh(sol));
            data.push_back(tmp_data);
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(3, 0.0, "sol_node", data, subdata);
        // Save the view
        nodedata.saveNodeData(filename, m_gmsh);
    }

    template<typename Function>
    void
    plot_tensor_function(const std::string& filename, const Function& fct) const
    {
        std::vector<gmsh::Data>    data;         // create data
        std::vector<gmsh::SubData> subdata;      // create subdata
        data.reserve(m_gmsh.getNumberofNodes()); // data has a size of nb_node

        auto storage = m_post_msh.mesh().backend_storage();

        // Compute the average value and save it
        for (int i_node = 1; i_node <= m_gmsh.getNumberofNodes(); i_node++)
        {
            const auto pt = storage->points[i_node - 1];

            const auto sol = fct(pt);

            const gmsh::Data tmp_data(i_node, disk::convertToVectorGmsh(sol));
            data.push_back(tmp_data);
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(9, 0.0, "sol_node", data, subdata);
        // Save the view
        nodedata.saveNodeData(filename, m_gmsh);
    }

    void
    save_vector_hho_solution_continuous(const std::string&              filename,
                                        const mesh_type&                msh,
                                        const hho_degree_info&          hdi,
                                        const std::vector<vector_type>& cells_solution) const
    {
        auto storage = m_post_msh.mesh().backend_storage();

        const static_vector<scalar_type, dimension> vzero = static_vector<scalar_type, dimension>::Zero();

        const auto cbs = vector_basis_size(hdi.cell_degree(), dimension, dimension);

        // first(number of data at this node), second(cumulated value)
        std::vector<std::pair<size_t, static_vector<scalar_type, dimension>>> value(m_gmsh.getNumberofNodes(),
                                                                                    std::make_pair(0, vzero));

        int cell_i = 0;
        for (auto& cl : msh)
        {
            const auto        cb         = make_vector_monomial_basis(msh, cl, hdi.cell_degree());
            const vector_type x          = cells_solution.at(cell_i).head(cbs);
            const auto        cell_nodes = m_post_msh.nodes_cell(cell_i);

            // Loop on the nodes of the cell
            for (auto& point_id : cell_nodes)
            {
                const auto pt = storage->points[point_id];

                const auto phi = cb.eval_functions(pt);
                const auto sol = eval(x, phi);

                // Add displacement at node
                value[point_id].first++;
                value[point_id].second += sol;
            }
            cell_i++;
        }

        std::vector<gmsh::Data>    data;         // create data
        std::vector<gmsh::SubData> subdata;      // create subdata
        data.reserve(m_gmsh.getNumberofNodes()); // data has a size of nb_node

        // Compute the average value and save it
        for (int i_node = 0; i_node < value.size(); i_node++)
        {
            const static_vector<scalar_type, dimension> sol_avr = value[i_node].second / double(value[i_node].first);

            const gmsh::Data tmp_data(i_node + 1, disk::convertToVectorGmsh(sol_avr));
            data.push_back(tmp_data);
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(3, 0.0, "sol_node_cont", data, subdata);
        // Save the view
        nodedata.saveNodeData(filename, m_gmsh);
    }

    void
    save_vector_hho_solution_discontinuous(const std::string&              filename,
                                           const mesh_type&                msh,
                                           const hho_degree_info&          hdi,
                                           const std::vector<vector_type>& cells_solution) const
    {
        gmsh::Gmesh gmsh(dimension);
        auto        storage = msh.backend_storage();

        std::vector<gmsh::Data>          data;    // create data (not used)
        const std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point

        const auto cbs = vector_basis_size(hdi.cell_degree(), dimension, dimension);

        int cell_i   = 0;
        int nb_nodes = 0;
        for (auto& cl : msh)
        {
            const auto                    cb         = make_vector_monomial_basis(msh, cl, hdi.cell_degree());
            const vector_type       x          = cells_solution.at(cell_i++).head(cbs);
            const auto              cell_nodes = disk::cell_nodes(msh, cl);
            std::vector<gmsh::Node> new_nodes;

            // loop on the nodes of the cell
            for (int i = 0; i < cell_nodes.size(); i++)
            {
                nb_nodes++;
                const auto point_ids = cell_nodes[i];
                const auto pt        = storage->points[point_ids];

                const auto phi = cb.eval_functions(pt);
                const auto sol = eval(x, phi);

                const std::vector<double>   sol_v = convertToVectorGmsh(sol);
                const std::array<double, 3> coor  = init_coor(pt);

                // Add a node
                const gmsh::Node tmp_node(coor, nb_nodes, 0);
                new_nodes.push_back(tmp_node);
                gmsh.addNode(tmp_node);

                const gmsh::Data datatmp(nb_nodes, sol_v);
                data.push_back(datatmp);
            }
            // Add new element
            disk::add_element(gmsh, new_nodes);
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(3, 0.0, "sol_node_disc", data, subdata);

        // Save the view
        nodedata.saveNodeData(filename, gmsh);
    }

    void
    save_vector_hho_deformed_continuous(const std::string&              filename,
                                           const mesh_type&                msh,
                                           const hho_degree_info&          hdi,
                                           const std::vector<vector_type>& cells_solution) const
    {
        gmsh::Gmesh gmsh(dimension);
        auto        storage = msh.backend_storage();

        const static_vector<scalar_type, dimension> vzero = static_vector<scalar_type, dimension>::Zero();
        const size_t                                nb_nodes(msh.points_size());
        const auto cbs = vector_basis_size(hdi.cell_degree(), dimension, dimension);

        // first(number of data at this node), second(cumulated value)
        std::vector<std::pair<size_t, static_vector<scalar_type, dimension>>> value(nb_nodes, std::make_pair(0, vzero));

        int cell_i = 0;
        for (auto& cl : msh)
        {
            const auto              cb         = disk::make_vector_monomial_basis(msh, cl, hdi.cell_degree());
            const vector_type x          = cells_solution.at(cell_i++).head(cbs);
            const auto        cell_nodes = disk::cell_nodes(msh, cl);

            // Loop on the nodes of the cell
            for (int i = 0; i < cell_nodes.size(); i++)
            {
                const auto point_ids = cell_nodes[i];
                const auto pt        = storage->points[point_ids];

                const auto phi  = cb.eval_functions(pt);
                const auto depl = eval(x, phi);

                // Add displacement at node
                value[point_ids].first++;
                value[point_ids].second += depl;
            }
        }

        // New coordinate
        int i_node = 0;
        for (auto itor = msh.points_begin(); itor != msh.points_end(); itor++)
        {
            const auto            pt   = *itor;
            std::array<double, 3> coor = init_coor(pt);

            const static_vector<scalar_type, dimension> depl_avr = value[i_node].second / double(value[i_node].first);

            for (int j = 0; j < dimension; j++)
                coor[j] += depl_avr(j);

            i_node++;
            const gmsh::Node tmp_node(coor, i_node, 0);
            gmsh.addNode(tmp_node);
        }
        const auto Nodes = gmsh.getNodes();

        // Add new elements
        for (auto& cl : msh)
        {
            const auto              cell_nodes = disk::cell_nodes(msh, cl);
            std::vector<gmsh::Node> new_nodes;

            // Loop on nodes of the cell
            for (int i = 0; i < cell_nodes.size(); i++)
            {
                const auto point_ids = cell_nodes[i];

                new_nodes.push_back(Nodes[point_ids]);
            }
            // Add new element
            disk::add_element(gmsh, new_nodes);
        }
        // Save mesh
        gmsh.writeGmesh(filename, 2);
    }

    void
    save_vector_hho_deformed_discontinuous(const std::string&              filename,
                                           const mesh_type&                msh,
                                           const hho_degree_info&          hdi,
                                           const std::vector<vector_type>& cells_solution) const
    {
        gmsh::Gmesh gmsh(dimension);
        auto        storage = msh.backend_storage();

        const auto cbs = vector_basis_size(hdi.cell_degree(), dimension, dimension);

        int    cell_i   = 0;
        size_t nb_nodes = 0;
        for (auto& cl : msh)
        {
            const auto                    cb         = make_vector_monomial_basis(msh, cl, hdi.cell_degree());
            const vector_type       x          = cells_solution.at(cell_i++).head(cbs);
            const auto              cell_nodes = disk::cell_nodes(msh, cl);
            std::vector<gmsh::Node> new_nodes;

            // Loop on nodes of the cell
            for (int i = 0; i < cell_nodes.size(); i++)
            {
                nb_nodes++;
                const auto point_ids = cell_nodes[i];
                const auto pt        = storage->points[point_ids];

                const auto phi  = cb.eval_functions(pt);
                const auto depl = eval(x, phi);

                std::array<double, 3> coor = disk::init_coor(pt);
                // Compute new coordinates
                for (int j = 0; j < dimension; j++)
                    coor[j] += depl(j);

                // Save node
                const gmsh::Node tmp_node(coor, nb_nodes, 0);
                new_nodes.push_back(tmp_node);
                gmsh.addNode(tmp_node);
            }
            // Add new element
            disk::add_element(gmsh, new_nodes);
        }
        // Save mesh
        gmsh.writeGmesh(filename, 2);
    }
};

} // end disk