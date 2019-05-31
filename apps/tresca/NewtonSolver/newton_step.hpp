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

// NewtonRaphson_step

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "../Informations.hpp"
#include "../Parameters.hpp"
#include "bases/bases.hpp"
#include "boundary_conditions/boundary_conditions.hpp"
#include "methods/hho"
#include "quadratures/quadratures.hpp"
#include "tresca_elementary_computation.hpp"

#include "solvers/solver.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

namespace MK
{

template<typename MeshType>
class NewtonRaphson_step_tresca
{
    typedef MeshType                                    mesh_type;
    typedef typename mesh_type::coordinate_type         scalar_type;
    typedef ParamRun<scalar_type>                       param_type;
    typedef typename disk::hho_degree_info              hdi_type;
    typedef disk::vector_boundary_conditions<mesh_type> bnd_type;
    typedef disk::MaterialData<scalar_type>             material_type;

    const static int dimension = mesh_type::dimension;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    typedef disk::assembler_mechanics<mesh_type> assembler_type;

    typedef tresca<mesh_type> elem_type;

    vector_type m_system_solution;

    const mesh_type&    m_msh;
    const hdi_type&     m_hdi;
    const bnd_type&     m_bnd;
    const param_type&   m_rp;
    const material_type m_material_data;
    assembler_type      m_assembler;

    std::vector<vector_type> m_bL;
    std::vector<matrix_type> m_AL;

    std::vector<vector_type> m_postprocess_data, m_solution, m_solution_faces;

    std::vector<bool> m_has_contact_faces;

    scalar_type m_F_int;

    bool m_verbose;

    void
    make_contact_face()
    {
        // cells with contact faces
        m_has_contact_faces.clear();
        m_has_contact_faces.assign(m_msh.cells_size(), false);

        size_t i = 0;
        for (auto& cl : m_msh)
        {
            m_has_contact_faces[i++] = m_bnd.cell_has_contact_faces(cl);
        }
    }

  public:
    NewtonRaphson_step_tresca(const mesh_type&    msh,
                              const hdi_type&     hdi,
                              const bnd_type&     bnd,
                              const param_type&   rp,
                              const material_type material_data) :
      m_msh(msh),
      m_hdi(hdi), m_rp(rp), m_bnd(bnd), m_verbose(rp.m_verbose), m_material_data(material_data)
    {
        m_AL.clear();
        m_AL.resize(m_msh.cells_size());

        m_bL.clear();
        m_bL.resize(m_msh.cells_size());

        m_assembler = disk::make_mechanics_assembler(m_msh, m_hdi, m_bnd);

        make_contact_face();
    }

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

    void
    initialize(const std::vector<vector_type>& initial_solution_faces, const std::vector<vector_type>& initial_solution)
    {
        m_solution_faces.clear();
        m_solution_faces = initial_solution_faces;
        assert(m_msh.faces_size() == m_solution_faces.size());

        m_solution.clear();
        m_solution = initial_solution;
        assert(m_msh.cells_size() == m_solution.size());
    }

    template<typename LoadFunction>
    AssemblyInfo
    assemble(const LoadFunction&             lf,
             const std::vector<matrix_type>& gradient_precomputed,
             const std::vector<matrix_type>& stab_precomputed)
    {
        elem_type    elem(m_msh, m_hdi, m_material_data, m_rp, m_bnd);
        AssemblyInfo ai;

        // set RHS to zero
        m_assembler.setZeroRhs();
        m_F_int = 0.0;

        timecounter tc, ttot;

        ttot.tic();
        int cell_i = 0;

        for (auto& cl : m_msh)
        {
            // Gradient Reconstruction
            matrix_type ET = gradient_precomputed[cell_i];

            // Stabilisation Contribution

            matrix_type stab;

            if (m_rp.m_stab)
            {
                stab = stab_precomputed.at(cell_i);
            }

            // Begin Assembly
            // Build rhs and lhs

            // Mechanical Computation

            tc.tic();
            elem.compute(cl, lf, ET, stab, m_solution.at(cell_i), m_has_contact_faces[cell_i]);

            matrix_type lhs = elem.K_int;
            vector_type rhs = elem.RTF;
            m_F_int += elem.F_int.squaredNorm();

            tc.toc();
            ai.m_time_elem += tc.to_double();
            ai.m_time_contact += elem.time_contact;

            // Static Condensation
            tc.tic();
            const auto scnp = make_vector_static_condensation_withMatrix(m_msh, cl, m_hdi, m_bnd, lhs, rhs);

            m_AL[cell_i] = std::get<1>(scnp);
            m_bL[cell_i] = std::get<2>(scnp);

            tc.toc();
            ai.m_time_statcond += tc.to_double();

            m_assembler.assemble_nl(m_msh, cl, m_bnd, std::get<0>(scnp), m_solution_faces);

            cell_i++;
        }

        m_F_int = sqrt(m_F_int);

        m_assembler.impose_neumann_boundary_conditions(m_msh, m_bnd);
        m_assembler.finalize();

        ttot.toc();
        ai.m_time_assembly      = ttot.to_double();
        ai.m_linear_system_size = m_assembler.LHS.rows();
        return ai;
    }

    SolveInfo
    solve(void)
    {
        timecounter tc;

        tc.tic();
        m_system_solution = vector_type::Zero(m_assembler.LHS.rows());

        disk::solvers::pardiso_params<scalar_type> pparams;
        mkl_pardiso(pparams, m_assembler.LHS, m_assembler.RHS, m_system_solution);
        tc.toc();

        return SolveInfo(m_assembler.LHS.rows(), m_assembler.LHS.nonZeros(), tc.to_double());
    }

    scalar_type
    postprocess()
    {
        timecounter tc;
        tc.tic();
        const auto cbs = disk::vector_basis_size(m_hdi.cell_degree(), dimension, dimension);

        const auto solF = m_assembler.expand_solution_nl(m_msh, m_bnd, m_system_solution, m_solution_faces);

        assert(m_solution_faces.size() == solF.size());
        // Update  unknowns
        // Update face Uf^{i+1} = Uf^i + delta Uf^i
        for (int i = 0; i < m_solution_faces.size(); i++)
        {
            // std::cout << i << " " << m_solution_faces.at(i).size() << " vs " << solF.at(i).size() << std::endl;

            assert(m_solution_faces.at(i).size() == solF.at(i).size());
            m_solution_faces.at(i) += solF.at(i);
        }

        // Update cell
        int cell_i = 0;
        for (auto& cl : m_msh)
        {
            const auto ci = disk::contact_info<mesh_type>(m_msh, cl, m_hdi, m_bnd);
            // Extract the solution
            const auto fcs            = ci.faces();
            const auto num_faces      = fcs.size();
            const auto num_faces_dofs = ci.num_faces_dofs();

            vector_type xFs    = vector_type::Zero(num_faces_dofs);
            size_t      offset = 0;
            for (int face_i = 0; face_i < num_faces; face_i++)
            {
                const auto face_id       = m_msh.lookup(fcs[face_i]);
                const auto solFi         = solF.at(face_id);
                const auto fbs           = solFi.size();
                xFs.segment(offset, fbs) = solFi;
                offset += fbs;
            }
            // Update element U^{i+1} = U^i + delta U^i ///
            m_solution.at(cell_i) +=
              disk::make_vector_static_decondensation_withMatrix(m_AL[cell_i], m_bL[cell_i], xFs);

            // std::cout << "KT_F " << m_AL[cell_i].norm() << std::endl;
            // std::cout << "sol_F" << std::endl;
            // std::cout << xFs.transpose() << std::endl;
            // std::cout << "ft" << std::endl;
            // std::cout << m_bL[cell_i].transpose() << std::endl;
            // std::cout << "sol_T" << std::endl;
            // std::cout << xT.transpose() << std::endl;
            // std::cout << (m_solution.at(cell_i)).segment(0, cbs).transpose() << std::endl;

            cell_i++;
        }

        tc.toc();
        return tc.to_double();
    }

    bool
    test_convergence(const int iter)
    {
        // norm of the solution
        scalar_type error_un = 0;
        for (int i = 0; i < m_solution_faces.size(); i++)
        {
            scalar_type norm = m_solution_faces[i].norm();
            error_un += norm * norm;
        }

        error_un = std::sqrt(error_un);

        if (error_un <= scalar_type(10E-15))
        {
            error_un = scalar_type(10E16);
        }

        // norm of the rhs
        const scalar_type residual  = m_assembler.RHS.norm();
        scalar_type       max_error = 0.0;
        for (int i = 0; i < m_assembler.RHS.size(); i++)
            max_error = std::max(max_error, std::abs(m_assembler.RHS(i)));

        // norm of the increment
        const scalar_type error_incr     = m_system_solution.norm();
        scalar_type       relative_displ = error_incr / error_un;
        scalar_type       relative_error = residual / m_F_int;

        if (iter == 0)
        {
            relative_displ = 1;
            relative_error = 1;
        }

        if (m_verbose)
        {
            std::string s_iter = "   " + std::to_string(iter) + "               ";
            s_iter.resize(9);

            if (iter == 0)
            {
                std::cout << "----------------------------------------------------------------------"
                             "------------------------"
                          << std::endl;
                std::cout << "| Iteration | Norme l2 incr | Relative incr |  Residual l2  | "
                             "Relative error | Maximum error |"
                          << std::endl;
                std::cout << "----------------------------------------------------------------------"
                             "------------------------"
                          << std::endl;
            }
            std::ios::fmtflags f(std::cout.flags());
            std::cout.precision(5);
            std::cout.setf(std::iostream::scientific, std::iostream::floatfield);
            std::cout << "| " << s_iter << " |   " << error_incr << " |   " << relative_displ << " |   " << residual
                      << " |   " << relative_error << "  |  " << max_error << "  |" << std::endl;
            std::cout << "-------------------------------------------------------------------------"
                         "---------------------"
                      << std::endl;
            std::cout.flags(f);
        }

        const scalar_type error = std::max(relative_displ, relative_error);

        if (error <= m_rp.m_epsilon)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    void
    save_solutions(std::vector<vector_type>& solution_faces, std::vector<vector_type>& solution)
    {
        solution_faces.clear();
        solution_faces = m_solution_faces;
        assert(m_solution_faces.size() == solution_faces.size());

        solution.clear();
        solution = m_solution;
        assert(m_solution.size() == solution.size());
    }
};

} // end NLE