/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2018                     nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

#pragma once

#include <vector>

#include "common/eigen.hpp"
#include "mechanics/behaviors/laws/law_cell_bones.hpp"
#include "mechanics/behaviors/maths_tensor.hpp"
#include "mechanics/behaviors/maths_utils.hpp"

namespace disk
{

// Law bones

template<typename MeshType, typename LawTypeQp, bool PlasticBehavior>
class LawTypeBones
{
  public:
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef LawTypeQp                           law_qp_type;
    typedef typename law_qp_type::data_type     data_type;

  private:
    typedef LawTypeCellBones<mesh_type, law_qp_type, PlasticBehavior> law_cell_type;

    size_t                     m_nb_qp;
    std::vector<law_cell_type> m_list_cell_qp;
    data_type                  m_data;

  public:
    LawTypeBones() : m_nb_qp(0){};

    LawTypeBones(const mesh_type& msh, const size_t degree)
    {
        m_nb_qp = 0;
        m_list_cell_qp.clear();
        m_list_cell_qp.reserve(msh.cells_size());

        for (auto& cl : msh)
        {
            law_cell_type cell_qp(msh, cl, degree);

            m_list_cell_qp.push_back(cell_qp);
            m_nb_qp += cell_qp.getNumberOfQP();
        }
    }

    void
    addMaterialData(const data_type materialData)
    {
        m_data = materialData;
    }

    data_type
    getMaterialData() const
    {
        return m_data;
    }

    int
    getNumberOfQP() const
    {
        return m_nb_qp;
    }

    void
    update()
    {
        for (auto& qp_cell : m_list_cell_qp)
        {
            qp_cell.update();
        }
    }

    law_cell_type&
    getCellQPs(const int cell_id)
    {
        return m_list_cell_qp.at(cell_id);
    }

    const law_cell_type&
    getCellQPs(const int cell_id) const
    {
        return m_list_cell_qp.at(cell_id);
    }

    law_cell_type
    getCellIVs(const int cell_id) const
    {
        return m_list_cell_qp.at(cell_id);
    }
};
}
