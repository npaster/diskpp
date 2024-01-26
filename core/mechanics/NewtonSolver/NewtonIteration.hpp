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
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

// NewtonRaphson iteration

#pragma once

#include <iostream>
#include <math.h>
#include <string>
#include <vector>

#include "NewtonSolverComput.hpp"
#include "NewtonSolverInformations.hpp"
#include "NewtonSolverParameters.hpp"
#include "TimeManager.hpp"
#include "bases/bases.hpp"

#include "adaptivity/adaptivity.hpp"
#include "boundary_conditions/boundary_conditions.hpp"
#include "mechanics/NewtonSolver/StabilizationManger.hpp"
#include "mechanics/behaviors/laws/behaviorlaws.hpp"
#include "mechanics/contact/ContactManager.hpp"

#include "methods/hho"

#include "solvers/solver.hpp"

#include "timecounter.h"

namespace disk
{

namespace mechanics
{

/**
 * @brief Newton-Raphson iteration for nonlinear solid mechanics
 *
 *  Specialized for HHO methods
 *
 *  Options :  - small and finite deformations
 *             - plasticity, hyperelasticity (various laws)
 *
 * @tparam MeshType type of the mesh
 */
template<typename MeshType>
class NewtonIteration
{
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    typedef NewtonSolverParameter<scalar_type>    param_type;
    typedef vector_boundary_conditions<mesh_type> bnd_type;
    typedef Behavior<mesh_type>                   behavior_type;

    typedef vector_mechanics_hho_assembler<mesh_type> assembler_type;
    typedef mechanical_computation<mesh_type>         elem_type;

    vector_type m_system_displ;

    assembler_type m_assembler;

    std::vector<vector_type> m_bL;
    std::vector<matrix_type> m_AL;

    std::vector<vector_type> m_postprocess_data;
    std::vector<vector_type> m_displ, m_displ_faces;
    std::vector<vector_type> m_velocity, m_acce;
    std::vector<vector_type> m_velocity_p, m_acce_p, m_acce_pred;

    TimeStep<scalar_type> m_time_step;

    scalar_type m_F_int;

    bool m_verbose;

  public:
    NewtonIteration(const mesh_type&                 msh,
                    const bnd_type&                  bnd,
                    const param_type&                rp,
                    const MeshDegreeInfo<mesh_type>& degree_infos,
                    const ContactManager<mesh_type>& contact_manager,
                    const TimeStep<scalar_type>&     current_step) :
      m_verbose(rp.m_verbose),
      m_time_step(current_step)
    {
        m_AL.clear();
        m_AL.resize(msh.cells_size());

        m_bL.clear();
        m_bL.resize(msh.cells_size());

        m_assembler = assembler_type(msh, degree_infos, bnd, contact_manager);
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
    initialize(const mesh_type&                msh,
               const param_type&               rp,
               const std::vector<vector_type>& initial_displ,
               const std::vector<vector_type>& initial_displ_faces,
               const std::vector<vector_type>& initial_velocity,
               const std::vector<vector_type>& initial_acce)
    {
        m_displ_faces.clear();
        m_displ_faces = initial_displ_faces;

        m_displ.clear();
        m_displ = initial_displ;

        m_velocity.clear();
        m_velocity = initial_velocity;

        m_acce.clear();
        m_acce = initial_acce;

        m_velocity_p.clear();
        m_velocity_p = initial_velocity;

        m_acce_p.clear();
        m_acce_p = initial_acce;

        if (rp.isUnsteady())
        {
            auto        dyna_rp = rp.getUnsteadyParameters();
            auto        beta    = dyna_rp["beta"];
            scalar_type dt      = m_time_step.increment_time();

            m_acce_pred.clear();
            m_acce_pred.reserve(msh.cells_size());

            for (auto& cl : msh)
            {
                const auto cl_id         = msh.lookup(cl);
                const auto num_cell_dofs = m_acce[cl_id].size();

                auto uT        = m_displ[cl_id].head(num_cell_dofs);
                auto acce_pred = -uT / (beta * dt * dt) - m_velocity[cl_id] / (beta * dt) -
                                 dt * dt / 2.0 * (1.0 - 2.0 * beta) * m_acce[cl_id];

                m_acce_pred.push_back(acce_pred);
            }
        }
    }

    template<typename LoadFunction>
    AssemblyInfo
    assemble(const mesh_type&                 msh,
             const bnd_type&                  bnd,
             const param_type&                rp,
             const MeshDegreeInfo<mesh_type>& degree_infos,
             const LoadFunction&              lf,
             const std::vector<matrix_type>&  gradient_precomputed,
             const std::vector<matrix_type>&  stab_precomputed,
             behavior_type&                   behavior,
             ContactManager<mesh_type>&       contact_manager,
             StabCoeffManager<scalar_type>&   stab_manager)
    {
        elem_type    elem;
        AssemblyInfo ai;

        // set RHS to zero
        m_assembler.initialize();
        m_F_int = 0.0;

        const bool small_def = (behavior.getDeformation() == SMALL_DEF);

        timecounter tc, ttot;

        ttot.tic();
        size_t cell_i = 0;

        for (auto& cl : msh)
        {
            const auto di = degree_infos.cellDegreeInfo(msh, cl);

            // Gradient Reconstruction
            // std::cout << "Grad" << std::endl;
            matrix_type GT;
            tc.tic();
            if (rp.m_precomputation)
            {
                GT = gradient_precomputed[cell_i];
            }
            else
            {
                if (small_def)
                {
                    const auto gradrec_sym_full = make_matrix_hho_symmetric_gradrec(msh, cl, degree_infos);
                    GT                          = gradrec_sym_full.first;
                }
                else
                {
                    const auto gradrec_full = make_matrix_hho_gradrec(msh, cl, degree_infos);
                    GT                      = gradrec_full.first;
                }
            }
            tc.toc();
            ai.m_time_gradrec += tc.to_double();

            // Begin Assembly
            // Build rhs and lhs

            // Mechanical Computation

            tc.tic();
            // std::cout << "Elem" << std::endl;
            elem.compute(msh,
                         cl,
                         rp,
                         degree_infos,
                         lf,
                         GT,
                         m_displ.at(cell_i),
                         m_acce_pred.at(cell_i),
                         m_time_step.increment_time(),
                         behavior,
                         stab_manager,
                         small_def);

            matrix_type lhs = elem.K_int;
            vector_type rhs = elem.RTF;
            m_F_int += elem.F_int.squaredNorm();

            tc.toc();
            ai.m_time_elem += tc.to_double();
            ai.m_time_law += elem.time_law;

            // Stabilisation Contribution
            // std::cout << "Stab" << std::endl;
            tc.tic();

            const auto beta = stab_manager.getValue(msh, cl);
            // std::cout << beta << std::endl;

            if (rp.m_stab)
            {
                if (rp.m_precomputation)
                {
                    const matrix_type& stab = stab_precomputed.at(cell_i);
                    assert(elem.K_int.rows() == stab.rows());
                    assert(elem.K_int.cols() == stab.cols());
                    assert(elem.RTF.rows() == stab.rows());
                    assert(elem.RTF.cols() == m_displ.at(cell_i).cols());

                    lhs += beta * stab;
                    rhs -= beta * stab * m_displ.at(cell_i);
                }
                else
                {
                    switch (rp.m_stab_type)
                    {
                        case HHO:
                        {
                            matrix_type stab_HHO;
                            // we do not make any difference for the displacement reconstruction
                            // if (small_def)
                            // {
                            //     const auto recons = make_vector_hho_symmetric_laplacian(msh, cl, degree_infos);
                            //     stab_HHO          = make_vector_hho_stabilization(msh, cl, recons.first,
                            //     degree_infos);
                            // }
                            // else
                            // {
                            const auto recons_scalar = make_scalar_hho_laplacian(msh, cl, degree_infos);
                            stab_HHO = make_vector_hho_stabilization_optim(msh, cl, recons_scalar.first, degree_infos);
                            // }

                            assert(elem.K_int.rows() == stab_HHO.rows());
                            assert(elem.K_int.cols() == stab_HHO.cols());
                            assert(elem.RTF.rows() == stab_HHO.rows());
                            assert(elem.RTF.cols() == m_displ.at(cell_i).cols());

                            lhs += beta * stab_HHO;
                            rhs -= beta * stab_HHO * m_displ.at(cell_i);
                            break;
                        }
                        case HDG:
                        {
                            const auto stab_HDG = make_vector_hdg_stabilization(msh, cl, degree_infos);
                            assert(elem.K_int.rows() == stab_HDG.rows());
                            assert(elem.K_int.cols() == stab_HDG.cols());
                            assert(elem.RTF.rows() == stab_HDG.rows());
                            assert(elem.RTF.cols() == m_displ.at(cell_i).cols());

                            lhs += beta * stab_HDG;
                            rhs -= beta * stab_HDG * m_displ.at(cell_i);
                            break;
                        }
                        case DG:
                        {
                            const auto stab_DG = make_vector_dg_stabilization(msh, cl, degree_infos);
                            assert(elem.K_int.rows() == stab_DG.rows());
                            assert(elem.K_int.cols() == stab_DG.cols());
                            assert(elem.RTF.rows() == stab_DG.rows());
                            assert(elem.RTF.cols() == m_displ.at(cell_i).cols());

                            lhs += beta * stab_DG;
                            rhs -= beta * stab_DG * m_displ.at(cell_i);
                            break;
                        }
                        case NO:
                        {
                            break;
                        }
                        default: throw std::invalid_argument("Unknown stabilization");
                    }
                }
            }
            tc.toc();
            ai.m_time_stab += tc.to_double();

            bool check_size = true;

            // // contact contribution
            // // std::cout << "Cont" << std::endl;
            // if (bnd.cell_has_contact_faces(cl))
            // {
            //     const auto cell_infos  = degree_infos.cellDegreeInfo(msh, cl);
            //     const auto faces_infos = cell_infos.facesDegreeInfo();

            //     const auto cell_degree = cell_infos.cell_degree();
            //     const auto grad_degree = cell_infos.grad_degree();

            //     const auto num_cell_dofs = vector_basis_size(cell_degree, mesh_type::dimension,
            //     mesh_type::dimension);

            //     const auto num_faces_dofs  = vector_faces_dofs(msh, faces_infos);
            //     const auto num_primal_dofs = num_cell_dofs + num_faces_dofs;
            //     const auto num_total_dofs  = num_primal_dofs;

            //     vector_type solution           = vector_type::Zero(num_total_dofs);
            //     solution.head(num_primal_dofs) = m_displ.at(cell_i);

            //     const auto fcs_cont = bnd.faces_with_contact(cl);

            //     size_t offset = num_primal_dofs;

            //     matrix_type Acont = matrix_type::Zero(num_total_dofs, num_total_dofs);
            //     vector_type rcont = vector_type::Zero(num_total_dofs);

            //     for (auto fc_cont : fcs_cont)
            //     {
            //         const auto fc_id   = msh.lookup(fc_cont);
            //         const auto mult_id = contact_manager.getMappingFaceToMult(fc_id);
            //         const auto face_id = contact_manager.getMappingMultToFace(mult_id);

            //         const auto num_mult_face = m_displ_mult.at(mult_id).size();

            //         solution.segment(offset, num_mult_face) = m_displ_mult.at(mult_id);
            //         rcont.segment(offset, num_mult_face)    = -m_displ_mult.at(mult_id);

            //         offset += num_mult_face;

            //         // std::cout << "id: " << fc_id << "->" << mult_id << "->" << face_id << std::endl;
            //         // std::cout << "mult: " << m_displ_mult.at(mult_id).transpose() << std::endl;
            //     }
            //     assert(offset == num_total_dofs);

            //     Acont.topLeftCorner(num_primal_dofs, num_primal_dofs) = lhs;
            //     rcont.head(num_primal_dofs)                           = rhs;

            //     Acont.bottomRightCorner(num_mult_dofs, num_mult_dofs) =
            //       matrix_type::Identity(num_mult_dofs, num_mult_dofs);

            //     // std::cout << "sol: " << solution.transpose() << std::endl;

            //     lhs = Acont;
            //     rhs = rcont;

            //     check_size = false;

            //     assert(lhs.rows() == num_total_dofs && lhs.cols() == num_total_dofs);
            //     assert(rhs.rows() == num_total_dofs);
            // }

            // Static Condensation
            // std::cout << "StatCond" << std::endl;
            tc.tic();
            const auto scnp = make_vector_static_condensation_withMatrix(msh, cl, degree_infos, lhs, rhs, check_size);

            m_AL[cell_i] = std::get<1>(scnp);
            m_bL[cell_i] = std::get<2>(scnp);

            tc.toc();
            ai.m_time_statcond += tc.to_double();

            const auto& lc = std::get<0>(scnp);
            // std::cout << "Assemb" << std::endl;
            m_assembler.assemble_nonlinear(msh, cl, bnd, contact_manager, lc.first, lc.second, m_displ_faces);

            cell_i++;
        }

        m_F_int = sqrt(m_F_int);
        // std::cout << "F_int: " << m_F_int << std::endl;

        m_assembler.impose_neumann_boundary_conditions(msh, bnd);
        m_assembler.finalize();

        ttot.toc();
        ai.m_time_assembly      = ttot.to_double();
        ai.m_linear_system_size = m_assembler.LHS.rows();
        return ai;
    }

    SolveInfo
    solve(void)
    {
        // std::cout << "begin solve" << std::endl;
        timecounter tc;

        tc.tic();
        m_system_displ = vector_type::Zero(m_assembler.LHS.rows());

        solvers::pardiso_params<scalar_type> pparams;
        mkl_pardiso(pparams, m_assembler.LHS, m_assembler.RHS, m_system_displ);
        tc.toc();

        // std::cout << "end solve" << std::endl;

        return SolveInfo(m_assembler.LHS.rows(), m_assembler.LHS.nonZeros(), tc.to_double());
    }

    scalar_type
    postprocess(const mesh_type&                 msh,
                const bnd_type&                  bnd,
                const param_type&                rp,
                const MeshDegreeInfo<mesh_type>& degree_infos,
                const ContactManager<mesh_type>& contact_manager)
    {
        // std::cout << "begin post_process" << std::endl;
        timecounter tc;
        tc.tic();

        // std::cout << m_system_displ << std::endl;

        // Update cell
        size_t cell_i = 0;
        for (auto& cl : msh)
        {
            const vector_type xdT =
              m_assembler.take_local_solution_nonlinear(msh, cl, bnd, m_system_displ, m_displ_faces);

            vector_type x_cond = xdT;

            // static decondensation
            const vector_type xT = m_bL[cell_i] - m_AL[cell_i] * x_cond;

            assert(m_displ.at(cell_i).size() == xT.size() + xdT.size());
            // Update element U^{i+1} = U^i + delta U^i ///
            m_displ.at(cell_i).head(xT.size()) += xT;
            m_displ.at(cell_i).segment(xT.size(), xdT.size()) += xdT;

            if (rp.isUnsteady())
            {
                auto        dyna_rp = rp.getUnsteadyParameters();
                auto        beta    = dyna_rp["beta"];
                auto        gamma   = dyna_rp["gamma"];
                scalar_type dt      = m_time_step.increment_time();

                m_acce.at(cell_i) = xT / (beta * dt * dt) + m_acce_pred.at(cell_i);
                m_velocity.at(cell_i) =
                  m_velocity_p.at(cell_i) + dt * ((1.0 - gamma) * m_acce.at(cell_i) + gamma * m_acce.at(cell_i));
            }

            // std::cout << "KT_F " << m_AL[cell_i].norm() << std::endl;
            // std::cout << "sol_F" << std::endl;
            // std::cout << xdT.transpose() << std::endl;
            // std::cout << "ft" << std::endl;
            // std::cout << m_bL[cell_i].transpose() << std::endl;
            // std::cout << "sol_T" << std::endl;
            // std::cout << xT.transpose() << std::endl;
            // std::cout << (m_displ.at(cell_i)).segment(0, xT.size()).transpose() << std::endl;

            cell_i++;
        }

        // Update  unknowns
        // Update face Uf^{i+1} = Uf^i + delta Uf^i
        size_t face_i = 0;
        for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
        {
            const auto fc = *itor;
            m_displ_faces.at(face_i++) +=
              m_assembler.take_local_solution_nonlinear(msh, fc, bnd, m_system_displ, m_displ_faces);
        }

        // std::cout << "end post_process" << std::endl;
        tc.toc();
        return tc.to_double();
    }

    bool
    convergence(const param_type& rp, const size_t iter)
    {
        // norm of the solution
        scalar_type error_un = 0;
        for (size_t i = 0; i < m_displ_faces.size(); i++)
        {
            scalar_type norm = m_displ_faces[i].norm();
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
        for (size_t i = 0; i < m_assembler.RHS.size(); i++)
            max_error = std::max(max_error, std::abs(m_assembler.RHS(i)));

        // norm of the increment
        const scalar_type error_incr     = m_system_displ.norm();
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

        if (!isfinite(error))
            throw std::runtime_error("Norm of residual is not finite");

        if (error <= rp.m_epsilon)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    void
    save_solutions(std::vector<vector_type>& displ,
                   std::vector<vector_type>& displ_faces,
                   std::vector<vector_type>& velocity,
                   std::vector<vector_type>& acce)
    {
        displ_faces.clear();
        displ_faces = m_displ_faces;
        assert(m_displ_faces.size() == displ_faces.size());

        displ.clear();
        displ = m_displ;
        assert(m_displ.size() == displ.size());

        velocity.clear();
        velocity = m_velocity;
        assert(m_velocity.size() == velocity.size());

        acce.clear();
        acce = m_acce;
        assert(m_acce.size() == acce.size());
    }
};
}
} // end diskpp