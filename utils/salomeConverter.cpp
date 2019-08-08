/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
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
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#include <algorithm>
#include <array>
#include <assert.h>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <unistd.h>
#include <vector>

class meshSalome
{

    class Edge
    {
      public:
        Edge() = delete;
        Edge(size_t ap0, size_t ap1, size_t iad)
        {
            if (ap0 < ap1)
            {
                p0 = ap0;
                p1 = ap1;
            }
            else if (ap0 > ap1)
            {
                p0 = ap1;
                p1 = ap0;
            }
            else
            {
                throw std::invalid_argument("an edge with the same vertices");
            }

            id = iad;
        }

        size_t              p0, p1;
        size_t              id;
        std::vector<size_t> faces;

        inline friend bool
        operator==(const Edge& lhs, const Edge& rhs)
        {
            if (lhs.p0 == rhs.p0 && lhs.p1 == rhs.p1)
            {
                return true;
            }

            return false;
        }
    };

    class Face
    {
      public:
        std::vector<size_t> vertices;
        std::vector<size_t> edges;
        size_t              id;
        int                 v0, v1;

        Face() = delete;
        Face(const std::vector<size_t>& vert, const size_t aid) : id(aid), v0(-1), v1(-1)
        {
            vertices.clear();
            vertices                          = vert;
            std::vector<size_t>::iterator min = std::min_element(std::begin(vertices), std::end(vertices));
            std::rotate(vertices.begin(), min, vertices.end());
        }

        void
        sort()
        {
            std::vector<size_t>::iterator min = std::min_element(std::begin(edges), std::end(edges));
            std::rotate(edges.begin(), min, edges.end());
        }

        inline friend bool
        operator==(const Face& lhs, const Face& rhs)
        {
            if (lhs.vertices.size() == rhs.vertices.size())
            {
                for (size_t i = 0; i < lhs.vertices.size(); i++)
                {
                    if (std::find(lhs.vertices.begin(), lhs.vertices.end(), rhs.vertices[i]) == lhs.vertices.end())
                    {
                        return false;
                    }
                }

                return true;
            }

            return false;
        }
    };

    class Volume
    {
      public:
        std::vector<size_t> vertices;
        std::vector<size_t> faces;
        size_t              id;

        Volume() = delete;
        Volume(const std::vector<size_t>& vert, const size_t iad)
        {
            vertices.clear();
            vertices                          = vert;
            std::vector<size_t>::iterator min = std::min_element(std::begin(vertices), std::end(vertices));
            std::rotate(vertices.begin(), min, vertices.end());
            id = iad;
        }

        void
        sort()
        {
            std::vector<size_t>::iterator min = std::min_element(std::begin(faces), std::end(faces));
            std::rotate(faces.begin(), min, faces.end());
        }

        inline friend bool
        operator==(const Volume& lhs, const Volume& rhs)
        {
            if (lhs.vertices.size() == rhs.vertices.size())
            {
                for (size_t i = 0; i < lhs.vertices.size(); i++)
                {
                    if (lhs.vertices[i] != rhs.vertices[i])
                    {
                        return false;
                    }
                }

                return true;
            }

            return false;
        }
    };

    std::vector<std::array<double, 3>> m_vertices;
    std::vector<Edge>                  m_edges;
    std::vector<Face>                  m_faces;
    std::vector<Volume>                m_volumes;
    std::vector<size_t>                m_bnd_faces;

    size_t ndim;

    size_t
    find(Edge& ed)
    {
        for (size_t i = 0; i < m_edges.size(); i++)
        {
            if (ed == m_edges[i])
            {
                m_edges[i].faces.push_back(ed.faces[0]);
                return i;
            }
        }

        m_edges.push_back(ed);

        return m_edges.size() - 1;
    }

    size_t
    find(Face& face)
    {
        for (size_t i = 0; i < m_faces.size(); i++)
        {
            if (face == m_faces[i])
            {
                if (m_faces[i].v0 == -1)
                {
                    m_faces[i].v0 = face.v0;
                }
                else if (m_faces[i].v1 == -1)
                {
                    m_faces[i].v1 = face.v0;
                }
                else
                {
                    throw std::invalid_argument("wrong face");
                }

                return i;
            }
        }

        // find edge
        size_t num_vert = face.vertices.size();
        size_t face_id  = m_faces.size();
        for (size_t i = 0; i < num_vert; i++)
        {
            Edge e0(face.vertices[i], face.vertices[(i + 1) % num_vert], 0);
            e0.faces.push_back(face_id);
            face.edges.push_back(find(e0));
        }

        m_faces.push_back(face);

        return face_id;
    }

    void
    search_boundaries()
    {
        m_bnd_faces.clear();

        if (ndim == 2)
        {
            for (size_t i = 0; i < m_edges.size(); i++)
            {
                if (m_edges[i].faces.size() == 1)
                {
                    m_bnd_faces.push_back(i);
                }
            }
        }
        else if (ndim == 3)
        {
            for (size_t i = 0; i < m_faces.size(); i++)
            {
                const auto fc = m_faces[i];

                if (fc.v0 == -1)
                {
                    throw std::invalid_argument("the face is not connected to cell");
                }

                if (fc.v1 == -1)
                {
                    m_bnd_faces.push_back(i);
                }
            }
        }
        else
        {
            throw std::invalid_argument("wrong dimension");
        }

        if (m_bnd_faces.empty())
        {
            throw std::invalid_argument("No boundary faces have been found. This is strange...");
        }

        std::cout << "Number of boundary faces: " << m_bnd_faces.size() << std::endl;
    }

    void
    check()
    {
        if (ndim == 2)
        {
            for (auto& vert : m_vertices)
            {
                if (std::abs(vert[2]) > 1E-13)
                {
                    throw std::invalid_argument("In 2D, the third coordinate of a vertice has to be zero");
                }
            }
        }

        for (Edge& e : m_edges)
        {
            if (e.p0 == e.p1)
            {
                throw std::invalid_argument("an edge with two identical vertices");
            }

            if (e.faces.size() == 0)
            {
                throw std::invalid_argument("an edge shared by no faces");
            }

            if (ndim == 2 && e.faces.size() > 2)
            {
                throw std::invalid_argument("an edge shared by more than two faces in 2D");
            }

            if (ndim == 3 && e.faces.size() == 1)
            {
                throw std::invalid_argument("an edge shared by only one face in 3D");
            }
        }

        for (Face& fc : m_faces)
        {
            for (size_t i = 0; i < fc.vertices.size(); i++)
                for (size_t j = i + 1; j < fc.vertices.size(); j++)
                    if (fc.vertices[i] == fc.vertices[j])
                    {
                        throw std::invalid_argument("a face with two identical vertices");
                    }

            for (size_t i = 0; i < fc.edges.size(); i++)
                for (size_t j = i + 1; j < fc.edges.size(); j++)
                    if (fc.edges[i] == fc.edges[j])
                    {
                        throw std::invalid_argument("a face with two identical edges");
                    }
        }

        for (Volume& vol : m_volumes)
        {
            for (size_t i = 0; i < vol.vertices.size(); i++)
                for (size_t j = i + 1; j < vol.vertices.size(); j++)
                    if (vol.vertices[i] == vol.vertices[j])
                    {
                        throw std::invalid_argument("a volume with two identical vertices");
                    }

            for (size_t i = 0; i < vol.faces.size(); i++)
                for (size_t j = i + 1; j < vol.faces.size(); j++)
                    if (vol.faces[i] == vol.faces[j])
                    {
                        throw std::invalid_argument("a volume with two identical faces");
                    }
        }
    }

  public:
    meshSalome() = default;

    bool
    readMeditmesh(const std::string& filename)
    {
        std::cout << "Guessed mesh format: MEDIT 3D" << std::endl;

        std::ifstream ifs(filename);
        std::string   keyword;

        if (!ifs.is_open())
        {
            std::cout << "Error opening " << filename << std::endl;
        }

        ifs >> keyword;
        if (keyword != "MeshVersionFormatted")
        {
            std::invalid_argument("Expected keyword \"MeshVersionFormatted\"");
        }

        size_t format;
        ifs >> format;

        if (format != 2)
        {
            std::cout << "Expected format 2 (here: " << format << ")" << std::endl;
            std::invalid_argument("Expected format");
        }

        ifs >> keyword;
        if (keyword != "Dimension")
        {
            std::invalid_argument("Expected keyword \"Dimension\"");
        }

        size_t dim;
        ifs >> dim;

        if (dim != 3)
        {
            std::cout << "Expected dimension = 3 (here: " << dim << ")" << std::endl;
            std::invalid_argument("Wrong dimension");
        }

        size_t elements_to_read = 0;
        ndim                    = 0;

        m_vertices.clear();
        m_edges.clear();
        m_faces.clear();
        m_volumes.clear();

        ifs >> keyword;
        while (keyword != "End")
        {
            std::cout << keyword << std::endl;
            if (keyword == "Vertices")
            {
                ifs >> elements_to_read;
                std::cout << elements_to_read << std::endl;

                m_vertices.reserve(elements_to_read);

                std::array<double, 3> vertice = {0.0, 0.0, 0.0};

                for (size_t i = 0; i < elements_to_read; i++)
                {
                    ifs >> vertice[0] >> vertice[1] >> vertice[2] >> keyword;
                    m_vertices.push_back(vertice);
                }
            }
            else if (keyword == "Edges")
            {
                ifs >> elements_to_read;
                std::cout << elements_to_read << std::endl;

                size_t p0, p1, id;

                ndim = std::max(ndim, size_t(1));

                m_edges.reserve(elements_to_read);

                for (size_t i = 0; i < elements_to_read; i++)
                {
                    ifs >> p0 >> p1 >> id;
                    Edge e(p0 - 1, p1 - 1, id);

                    m_edges.push_back(e);
                }
            }
            else if (keyword == "Triangles")
            {
                ifs >> elements_to_read;
                std::cout << elements_to_read << std::endl;

                m_faces.reserve(m_faces.size() + elements_to_read);

                size_t p0, p1, p2, id;

                ndim = std::max(ndim, size_t(2));

                for (size_t i = 0; i < elements_to_read; i++)
                {
                    ifs >> p0 >> p1 >> p2 >> id;
                    Face tri(std::vector<size_t>({p0 - 1, p1 - 1, p2 - 1}), id);
                    find(tri);
                }
            }
            else if (keyword == "Quadrilaterals")
            {
                ifs >> elements_to_read;
                std::cout << elements_to_read << std::endl;

                m_faces.reserve(m_faces.size() + elements_to_read);

                size_t p0, p1, p2, p3, id;

                ndim = std::max(ndim, size_t(2));

                for (size_t i = 0; i < elements_to_read; i++)
                {
                    ifs >> p0 >> p1 >> p2 >> p3 >> id;
                    Face q(std::vector<size_t>({p0 - 1, p1 - 1, p2 - 1, p3 - 1}), id);
                    find(q);
                }
            }
            else if (keyword == "Tetrahedra")
            {
                ifs >> elements_to_read;
                std::cout << elements_to_read << std::endl;

                m_volumes.reserve(m_volumes.size() + elements_to_read);
                m_faces.reserve(m_faces.size() + 3 * elements_to_read);

                ndim = std::max(ndim, size_t(3));

                size_t p0, p1, p2, p3, id;
                for (size_t i = 0; i < elements_to_read; i++)
                {
                    ifs >> p0 >> p1 >> p2 >> p3 >> id;
                    Volume tet(std::vector<size_t>({p0 - 1, p1 - 1, p2 - 1, p3 - 1}), id);
                    size_t tet_id = m_volumes.size();

                    Face tri1(std::vector<size_t>({p0 - 1, p1 - 1, p2 - 1}), 0);
                    tri1.v0 = tet_id;
                    tet.faces.push_back(find(tri1));

                    Face tri2(std::vector<size_t>({p0 - 1, p3 - 1, p2 - 1}), 0);
                    tri2.v0 = tet_id;
                    tet.faces.push_back(find(tri2));

                    Face tri3(std::vector<size_t>({p1 - 1, p2 - 1, p3 - 1}), 0);
                    tri3.v0 = tet_id;
                    tet.faces.push_back(find(tri3));

                    Face tri4(std::vector<size_t>({p0 - 1, p3 - 1, p1 - 1}), 0);
                    tri4.v0 = tet_id;
                    tet.faces.push_back(find(tri4));

                    m_volumes.push_back(tet);
                }
            }
            else if (keyword == "Hexahedra")
            {
                ifs >> elements_to_read;
                std::cout << elements_to_read << std::endl;

                m_volumes.reserve(m_volumes.size() + elements_to_read);
                m_faces.reserve(m_faces.size() + 4 * elements_to_read);

                ndim = std::max(ndim, size_t(3));

                size_t p0, p1, p2, p3, p4, p5, p6, p7, id;
                for (size_t i = 0; i < elements_to_read; i++)
                {
                    ifs >> p0 >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> id;
                    Volume hex(std::vector<size_t>({p0 - 1, p1 - 1, p2 - 1, p3 - 1, p4 - 1, p5 - 1, p6 - 1, p7 - 1}),
                               id);
                    size_t hexid = m_volumes.size();

                    Face q1(std::vector<size_t>({p0 - 1, p1 - 1, p2 - 1, p3 - 1}), 0);
                    q1.v0 = hexid;
                    hex.faces.push_back(find(q1));

                    Face q2(std::vector<size_t>({p0 - 1, p4 - 1, p7 - 1, p3 - 1}), 0);
                    q2.v0 = hexid;
                    hex.faces.push_back(find(q2));

                    Face q3(std::vector<size_t>({p4 - 1, p5 - 1, p6 - 1, p7 - 1}), 0);
                    q3.v0 = hexid;
                    hex.faces.push_back(find(q3));

                    Face q4(std::vector<size_t>({p5 - 1, p1 - 1, p2 - 1, p6 - 1}), 0);
                    q4.v0 = hexid;
                    hex.faces.push_back(find(q4));

                    Face q5(std::vector<size_t>({p3 - 1, p7 - 1, p6 - 1, p2 - 1}), 0);
                    q5.v0 = hexid;
                    hex.faces.push_back(find(q5));

                    Face q6(std::vector<size_t>({p0 - 1, p4 - 1, p5 - 1, p1 - 1}), 0);
                    q6.v0 = hexid;
                    hex.faces.push_back(find(q6));

                    m_volumes.push_back(hex);
                }
            }
            else if (keyword == "Prisms")
            {
                ifs >> elements_to_read;
                std::cout << elements_to_read << std::endl;

                m_volumes.reserve(m_volumes.size() + elements_to_read);
                m_faces.reserve(m_faces.size() + 3 * elements_to_read);

                ndim = std::max(ndim, size_t(3));

                size_t p0, p1, p2, p3, p4, p5, id;
                for (size_t i = 0; i < elements_to_read; i++)
                {
                    ifs >> p0 >> p1 >> p2 >> p3 >> p4 >> p5 >> id;
                    Volume pri(std::vector<size_t>({p0 - 1, p1 - 1, p2 - 1, p3 - 1, p4 - 1, p5 - 1}), id);
                    size_t priid = m_volumes.size();

                    Face t1(std::vector<size_t>({p0 - 1, p1 - 1, p2 - 1}), 0);
                    t1.v0 = priid;
                    pri.faces.push_back(find(t1));

                    Face q1(std::vector<size_t>({p0 - 1, p3 - 1, p4 - 1, p1 - 1}), 0);
                    q1.v0 = priid;
                    pri.faces.push_back(find(q1));

                    Face q2(std::vector<size_t>({p2 - 1, p1 - 1, p4 - 1, p5 - 1}), 0);
                    q2.v0 = priid;
                    pri.faces.push_back(find(q2));

                    Face q3(std::vector<size_t>({p0 - 1, p2 - 1, p5 - 1, p3 - 1}), 0);
                    q3.v0 = priid;
                    pri.faces.push_back(find(q3));

                    Face t2(std::vector<size_t>({p3 - 1, p4 - 1, p5 - 1}), 0);
                    t2.v0 = priid;
                    pri.faces.push_back(find(t2));

                    m_volumes.push_back(pri);
                }
            }
            else
            {
                std::invalid_argument("Error parsing Medit file");
            }

            ifs >> keyword;
        }
        ifs.close();

        assert(ndim > 1);

        m_vertices.shrink_to_fit();
        m_edges.shrink_to_fit();
        m_faces.shrink_to_fit();
        m_volumes.shrink_to_fit();

        // order element

        for (auto& fc : m_faces)
        {
            fc.sort();
        }

        for (auto& vol : m_volumes)
        {
            vol.sort();
        }

        check();

        search_boundaries();

        std::cout << "The mesh is " << ndim << "D." << std::endl;

        return true;
    }

    void
    writeMeditmesh(const std::string& filename)
    {
        std::ofstream mesh_file(filename);
        if (!mesh_file.is_open())
        {
            throw std::invalid_argument("Unable to open file ");
        }

        mesh_file << "MeshVersionFormatted 2" << std::endl;
        mesh_file << "Dimension"
                  << " " << ndim << std::endl;
        mesh_file << " " << std::endl;

        mesh_file << "Vertices"
                  << "\t" << m_vertices.size() << std::endl;

        for (auto& vert : m_vertices)
        {
            if (ndim == 3)
            {
                mesh_file << vert[0] << " " << vert[1] << " " << vert[2] << std::endl;
            }
            else if (ndim == 2)
            {
                mesh_file << vert[0] << " " << vert[1] << std::endl;
            }
            else
            {
                throw std::invalid_argument("wrong dimension");
            }
        }

        if (ndim == 3)
        {
            mesh_file << "Volumes->Faces"
                      << "\t" << m_volumes.size() << std::endl;

            for (auto& vol : m_volumes)
            {
                const auto faces = vol.faces;
                mesh_file << faces.size();

                for (size_t i = 0; i < faces.size(); i++)
                {
                    mesh_file << " " << faces[i];
                }
                mesh_file << std::endl;
            }

            mesh_file << "Volumes->Vertices"
                      << "\t" << m_volumes.size() << std::endl;

            for (auto& vol : m_volumes)
            {
                mesh_file << vol.vertices.size();
                for (auto& fc : vol.vertices)
                {
                    mesh_file << " " << fc;
                }
                mesh_file << std::endl;
            }
        }

        mesh_file << "Faces->Edges"
                  << "\t" << m_faces.size() << std::endl;

        for (auto& fc : m_faces)
        {
            mesh_file << fc.edges.size();
            for (auto& e : fc.edges)
            {
                mesh_file << " " << e;
            }
            mesh_file << std::endl;
        }

        mesh_file << "Faces->Vertices"
                  << "\t" << m_faces.size() << std::endl;

        for (auto& fc : m_faces)
        {
            mesh_file << fc.vertices.size();
            for (auto& e : fc.vertices)
            {
                mesh_file << " " << e;
            }
            mesh_file << std::endl;
        }

        mesh_file << "Faces->Control volumes"
                  << "\t" << m_bnd_faces.size() << std::endl;

        for (size_t i = 0; i < m_bnd_faces.size(); i++)
        {
            mesh_file << m_bnd_faces[i] << " " << m_faces[m_bnd_faces[i]].id << std::endl;
        }

        mesh_file << "Edges"
                  << "\t" << m_edges.size() << std::endl;

        for (auto& edge : m_edges)
        {
            mesh_file << edge.p0 << " " << edge.p1 << std::endl;
        }

        mesh_file << "End" << std::endl;

        mesh_file.close();
    }

    void
    writePolymesh(const std::string& filename)
    {
        std::ofstream mesh_file(filename);
        if (!mesh_file.is_open())
        {
            throw std::invalid_argument("Unable to open file ");
        }

        mesh_file << "MeshVersionFormatted 2" << std::endl;
        mesh_file << "Dimension"
                  << " " << ndim << std::endl;
        mesh_file << " " << std::endl;

        mesh_file << "Vertices"
                  << "\t" << m_vertices.size() << std::endl;

        for (auto& vert : m_vertices)
        {
            if (ndim == 3)
            {
                mesh_file << vert[0] << " " << vert[1] << " " << vert[2] << std::endl;
            }
            else if (ndim == 2)
            {
                mesh_file << vert[0] << " " << vert[1] << std::endl;
            }
            else
            {
                throw std::invalid_argument("wrong dimension");
            }
        }

        mesh_file << "Edges"
                  << "\t" << m_edges.size() << std::endl;

        for (auto& edge : m_edges)
        {
            mesh_file << edge.p0 << " " << edge.p1 << std::endl;
        }

        mesh_file << "Faces->Edges"
                  << "\t" << m_faces.size() << std::endl;

        for (auto& fc : m_faces)
        {
            mesh_file << fc.edges.size();
            for (auto& e : fc.edges)
            {
                mesh_file << " " << e;
            }
            mesh_file << std::endl;
        }

        mesh_file << "Faces->Vertices"
                  << "\t" << m_faces.size() << std::endl;

        for (auto& fc : m_faces)
        {
            mesh_file << fc.vertices.size();
            for (auto& e : fc.vertices)
            {
                mesh_file << " " << e;
            }
            mesh_file << std::endl;
        }

        if (ndim == 3)
        {
            mesh_file << "Volumes->Faces"
                      << "\t" << m_volumes.size() << std::endl;

            for (auto& vol : m_volumes)
            {
                const auto faces = vol.faces;
                mesh_file << faces.size();

                for (size_t i = 0; i < faces.size(); i++)
                {
                    mesh_file << " " << faces[i];
                }
                mesh_file << std::endl;
            }

            mesh_file << "Volumes->Vertices"
                      << "\t" << m_volumes.size() << std::endl;

            for (auto& vol : m_volumes)
            {
                mesh_file << vol.vertices.size();
                for (auto& fc : vol.vertices)
                {
                    mesh_file << " " << fc;
                }
                mesh_file << std::endl;
            }
        }

        mesh_file << "Faces->Control volumes"
                  << "\t" << m_bnd_faces.size() << std::endl;

        for (size_t i = 0; i < m_bnd_faces.size(); i++)
        {
            mesh_file << m_bnd_faces[i] << " " << m_faces[m_bnd_faces[i]].id << std::endl;
        }

        mesh_file << "End" << std::endl;

        mesh_file.close();
    }

    void
    writeMsh(const std::string& filename)
    {
        std::ofstream mesh_file(filename);
        if (!mesh_file.is_open())
        {
            std::cout << "Unable to open file " << filename << std::endl;
            abort();
        }

        mesh_file << "MeshVersionFormatted 2" << std::endl;
        mesh_file << " " << std::endl;
        mesh_file << "Dimension 3" << std::endl;
        mesh_file << " " << std::endl;

        mesh_file << "Vertices" << std::endl;
        mesh_file << m_vertices.size() << std::endl;

        for (auto& vert : m_vertices)
        {
            mesh_file << vert[0] << " " << vert[1] << " " << vert[2] << " " << 0 << std::endl;
        }
        mesh_file << " " << std::endl;

        size_t n_tri = 0, n_quad = 0;

        mesh_file << "Edges" << std::endl;
        mesh_file << m_edges.size() << std::endl;

        for (auto& edge : m_edges)
        {
            mesh_file << edge.p0 + 1 << " " << edge.p1 + 1 << " " << edge.id << std::endl;
        }
        mesh_file << " " << std::endl;

        for (auto& fc : m_faces)
        {
            switch (fc.vertices.size())
            {
                case 3: n_tri++; break;

                case 4: n_quad++; break;

                default: throw std::invalid_argument("wrong number of vertices");
            }
        }

        if (n_tri > 0)
        {
            mesh_file << "Triangles" << std::endl;
            mesh_file << n_tri << std::endl;

            for (auto& fc : m_faces)
            {
                const auto vert = fc.vertices;

                if (vert.size() == 3)
                {
                    mesh_file << vert[0] + 1 << " " << vert[1] + 1 << " " << vert[2] + 1 << " " << fc.id << std::endl;
                    n_tri--;
                }
                assert(n_tri == 0);
            }
            mesh_file << " " << std::endl;
        }

        if (n_quad > 0)
        {
            mesh_file << "Quadrilaterals" << std::endl;
            mesh_file << n_quad << std::endl;

            for (auto& fc : m_faces)
            {
                const auto vert = fc.vertices;

                if (vert.size() == 4)
                {
                    mesh_file << vert[0] + 1 << " " << vert[1] + 1 << " " << vert[2] + 1 << " " << vert[3] + 1 << " "
                              << fc.id << std::endl;
                    n_quad--;
                }
            }
            assert(n_quad == 0);
            mesh_file << " " << std::endl;
        }

        size_t n_te = 0, n_prism = 0, n_hexa = 0;

        for (auto& vol : m_volumes)
        {
            switch (vol.vertices.size())
            {
                case 4: n_te++; break;

                case 6: n_prism++; break;

                case 8: n_hexa++; break;

                default: throw std::invalid_argument("wrong number of vertices");
            }
        }

        if (n_te > 0)
        {
            mesh_file << "Tetrahedra" << std::endl;
            mesh_file << n_te << std::endl;

            for (auto& vol : m_volumes)
            {
                const auto vert = vol.vertices;

                if (vert.size() == 4)
                {
                    mesh_file << vert[0] + 1 << " " << vert[1] + 1 << " " << vert[2] + 1 << " " << vert[3] + 1 << " "
                              << vol.id << std::endl;
                    n_te--;
                }
            }
            assert(n_te == 0);
            mesh_file << " " << std::endl;
        }

        if (n_hexa > 0)
        {
            mesh_file << "Hexahedra" << std::endl;
            mesh_file << n_hexa << std::endl;

            for (auto& vol : m_volumes)
            {
                const auto vert = vol.vertices;

                if (vert.size() == 8)
                {
                    mesh_file << vert[0] + 1 << " " << vert[1] + 1 << " " << vert[2] + 1 << " " << vert[3] + 1 << " "
                              << vert[4] + 1 << " " << vert[5] + 1 << " " << vert[6] + 1 << " " << vert[7] + 1 << " "
                              << vol.id << std::endl;
                    n_hexa--;
                }
            }
            assert(n_hexa == 0);
            mesh_file << " " << std::endl;
        }

        if (n_prism > 0)
        {
            mesh_file << "Prisms" << std::endl;
            mesh_file << n_prism << std::endl;

            for (auto& vol : m_volumes)
            {
                const auto vert = vol.vertices;

                if (vert.size() == 6)
                {
                    mesh_file << vert[0] + 1 << " " << vert[1] + 1 << " " << vert[2] + 1 << " " << vert[3] + 1 << " "
                              << vert[4] + 1 << " " << vert[5] + 1 << " " << vol.id << std::endl;
                    n_prism--;
                }
            }
            assert(n_prism == 0);
            mesh_file << " " << std::endl;
        }

        mesh_file << "End" << std::endl;

        mesh_file.close();
    }
};

/* To use:
 * ./utils/createnonconforming.cpp -n number_hexagons input_mesh.*  output_mesh.medit2d
 */

int
main(int argc, char** argv)
{
    argc -= optind;
    argv += optind;

    char* filename = argv[0];
    char* output   = argv[1];

    if (std::regex_match(filename, std::regex(".*\\.mesh$")))
    {
        meshSalome mesh;

        mesh.readMeditmesh(filename);
        mesh.writeMeditmesh(output);
        mesh.writeMsh("test_edge.mesh");
    }
    else
    {
        std::cout << "Unknown format" << std::endl;
    }
}

// lire maillage
// modifier maillage
// ecrire maillage
