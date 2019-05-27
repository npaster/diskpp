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
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

#pragma once

#include <vector>

#include "common/eigen.hpp"
#include "core/mechanics/behaviors/logarithmic_strain/LogarithmicStrain_qp.hpp"
#include "mechanics/behaviors/laws/behaviorlaws.hpp"
#include "mechanics/behaviors/laws/law_cell_bones.hpp"

namespace disk
{

namespace mechanics
{

// Routine for Logarithmic Stain

/* For details see the paper:
 *   Anisotropic additive plasticity in the logaritmic strain: modular
 *   kinematic formulation and implementation based on incremental minimization
 *   principles for standard materials
 *   C. Miehe, N. Apel, M. Lambrecht
 *   Comput. Methods Appl. Mech. Engrg. (2002)
 */

template<typename LawType>
class LogarithmicStrain
{
  public:
    typedef LawType                             law_hpp_type;
    typedef typename law_hpp_type::law_qp_type  law_hpp_qp_type;
    typedef typename law_hpp_type::mesh_type    mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename law_hpp_qp_type::data_type data_type;

    typedef LogarithmicStrain_qp<law_hpp_qp_type> law_qp_type;

    typedef LawTypeCellBones<mesh_type, law_qp_type, true> law_cell_type;

  private:
    size_t                     m_nb_qp;
    std::vector<law_cell_type> m_list_cell_qp;
    data_type                  m_data;

  public:
    LogarithmicStrain() : m_nb_qp(0){};

    LogarithmicStrain(const mesh_type& msh, const size_t degree)
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
    addMaterialData(const data_type& material_data)
    {
        m_data = material_data;
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

    law_cell_type
    getCellIVs(const int cell_id) const
    {
        return m_list_cell_qp.at(cell_id);
    }

    scalar_type
    smallest_eigenvalue() const
    {
        using std::min;

        scalar_type ev = 10.0E50;

        for (auto& qp_cell : m_list_cell_qp)
        {
            const auto& law_quadpoints = qp_cell.getQPs();

            for(auto& qp : law_quadpoints)
            {
                ev = min(ev, qp.small_eigenvalue());
            }
        }
        return ev;
    }

    void
    save_smallest_eigenvalues(const std::string& filename) const
    {
        std::ofstream output;
        output.open(filename, std::ofstream::out | std::ofstream::trunc);

        if (!output.is_open())
        {
            std::cerr << "Unable to open file " << filename << std::endl;
        }

        for (auto& qp_cell : m_list_cell_qp)
        {
            const auto& law_quadpoints = qp_cell.getQPs();

            for(auto& qp : law_quadpoints)
            {
                output << qp.smallest_eigenvalue() << std::endl;
            }
        }

        output.close();
    }

    void
    stat_smallest_eigenvalues() const
    {
        std::vector<scalar_type> tab;
        tab.reserve(m_nb_qp);

        for (auto& qp_cell : m_list_cell_qp)
        {
            const auto& law_quadpoints = qp_cell.getQPs();

            for(auto& qp : law_quadpoints)
            {
                tab.push_back(qp.small_eigenvalue());
            }
        }

        const auto small_ev = smallest_eigenvalue();

        int num_items;

        std::cout << "Statistique about eigenvalues: " << std::endl;
        num_items = std::count_if(tab.begin(), tab.end(), [](scalar_type ev){return ev <= -10.0;});
        std::cout << "** Eigenvalues between min < ev < -10: " << num_items << "  ( " << scalar_type(num_items)/m_nb_qp << " % )" << std::endl;
        num_items = std::count_if(tab.begin(), tab.end(), [](scalar_type ev){ if(ev > -10 && ev <= -5.0){return true;} return false;});
        std::cout << "** Eigenvalues between -10 < ev < -5 : " << num_items << "  ( " << scalar_type(num_items)/m_nb_qp << " % )" << std::endl;
        num_items = std::count_if(tab.begin(), tab.end(), [](scalar_type ev){ if(ev > -5 && ev <= -1.0){return true;} return false;});
        std::cout << "** Eigenvalues between -5 < ev < -1  : " << num_items << "  ( " << scalar_type(num_items)/m_nb_qp << " % )" << std::endl;
        num_items = std::count_if(tab.begin(), tab.end(), [](scalar_type ev){ if(ev > -1.0 && ev <= -0.5){return true;} return false;});
        std::cout << "** Eigenvalues between -1 < ev < -0.5: " << num_items << "  ( " << scalar_type(num_items)/m_nb_qp << " % )" << std::endl;
        num_items = std::count_if(tab.begin(), tab.end(), [](scalar_type ev){ if(ev > -0.5 && ev <= 0.0){return true;} return false;});
        std::cout << "** Eigenvalues between -0.5 < ev < 0 : " << num_items << "  ( " << scalar_type(num_items)/m_nb_qp << " % )" << std::endl;
        num_items = std::count_if(tab.begin(), tab.end(), [](scalar_type ev){ if(ev > 0.0 && ev <= 0.5){return true;} return false;});
        std::cout << "** Eigenvalues between 0 < ev < 0.5  : " << num_items << "  ( " << scalar_type(num_items)/m_nb_qp << " % )" << std::endl;
        num_items = std::count_if(tab.begin(), tab.end(), [](scalar_type ev){ if(ev > .5 && ev <= 1.0){return true;} return false;});
        std::cout << "** Eigenvalues between 0.5 < ev < 1  : " << num_items << "  ( " << scalar_type(num_items)/m_nb_qp << " % )" << std::endl;
        num_items = std::count_if(tab.begin(), tab.end(), [](scalar_type ev){ if(ev > 1.0 && ev <= 5.0){return true;} return false;});
        std::cout << "** Eigenvalues between 1 < ev < 5    : " << num_items << "  ( " << scalar_type(num_items)/m_nb_qp << " % )" << std::endl;
        num_items = std::count_if(tab.begin(), tab.end(), [](scalar_type ev){ if(ev > 5 && ev <= 10.0){return true;} return false;});
        std::cout << "** Eigenvalues between 5 < ev < 10   : " << num_items << "  ( " << scalar_type(num_items)/m_nb_qp << " % )" << std::endl;
        num_items = std::count_if(tab.begin(), tab.end(), [](scalar_type ev){ if(ev > 10 && ev <= 50.0){return true;} return false;});
        std::cout << "** Eigenvalues between 10 < ev < 50  : " << num_items << "  ( " << scalar_type(num_items)/m_nb_qp << " % )" << std::endl;
        num_items = std::count_if(tab.begin(), tab.end(), [](scalar_type ev){ if(ev > 50 ){return true;} return false;});
        std::cout << "** Eigenvalues between 50 < ev       : " << num_items << "  ( " << scalar_type(num_items)/m_nb_qp << " % )" << std::endl;

    }
};
}
}