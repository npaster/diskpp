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

#include <string>

namespace disk
{

template<typename mesh_type>
class mesh_loader
{
    bool m_verbose;

  public:
    mesh_loader() : m_verbose(false) {}

    virtual bool
    read_mesh(const std::string&)
    {
        return false;
    }
    virtual bool populate_mesh(mesh_type&) = 0;

    bool
    verbose(void) const
    {
        return m_verbose;
    }
    void
    verbose(bool v)
    {
        m_verbose = v;
    }

    virtual ~mesh_loader() {}
};

}