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

#include "diskpp/adaptivity/adaptivity.hpp"
#include "diskpp/bases/bases.hpp"
#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/common/timecounter.hpp"

#include "diskpp/mechanics/NewtonSolver/Fields.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolverComput.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolverDynamic.hpp"
#include "diskpp/mechanics/NewtonSolver/NewtonSolverInformations.hpp"
#include "diskpp/mechanics/NewtonSolver/NonLinearParameters.hpp"
#include "diskpp/mechanics/NewtonSolver/StabilizationManager.hpp"
#include "diskpp/mechanics/NewtonSolver/TimeManager.hpp"
#include "diskpp/mechanics/behaviors/laws/behaviorlaws.hpp"

#include "diskpp/methods/hho"
#include "diskpp/solvers/solver.hpp"

namespace disk {

namespace mechanics {

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
template <typename MeshType> class SplittingIteration {
    typedef MeshType mesh_type;
    typedef typename mesh_type::cell cell_type;
    typedef typename mesh_type::coordinate_type scalar_type;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    typedef NonLinearParameters<scalar_type> param_type;
    typedef vector_boundary_conditions<mesh_type> bnd_type;
    typedef Behavior<mesh_type> behavior_type;

    typedef vector_mechanics_hho_assembler<mesh_type> assembler_type;
    typedef mechanical_computation<mesh_type> elem_type;
    typedef dynamic_computation<mesh_type> dyna_type;

    typedef Eigen::PardisoLDLT<Eigen::SparseMatrix<scalar_type>> solver_type;

    vector_type m_system_displ, m_system_displ_prev, m_residual;

    assembler_type m_assembler;

    std::vector<vector_type> m_stab_uT;
    std::vector<matrix_type> m_lhs, m_stab_zFF;

    TimeStep<scalar_type> m_time_step;

    dyna_type m_dyna;

    bool m_verbose;

    solver_type m_solver;

    matrix_type _mass_term(const mesh_type &msh, const cell_type &cl,
                           const MeshDegreeInfo<mesh_type> &degree_infos) const {

        const auto cell_infos = degree_infos.cellDegreeInfo(msh, cl);
        const auto faces_infos = cell_infos.facesDegreeInfo();
        const auto num_faces_dofs = vector_faces_dofs(msh, faces_infos);

        matrix_type mm = matrix_type::Zero(num_faces_dofs, num_faces_dofs);

        const auto fcs = faces(msh, cl);

        auto to_vector = [](const matrix_type &scalar_matrix) {
            const int scal_total_dofs = scalar_matrix.rows();
            const int vect_total_tofs = scal_total_dofs * mesh_type::dimension;

            matrix_type mm = matrix_type::Zero(vect_total_tofs, vect_total_tofs);

            for (int i = 0; i < scal_total_dofs; i++) {
                const auto row = i * mesh_type::dimension;
                for (int j = 0; j < scal_total_dofs; j++) {
                    const auto col = j * mesh_type::dimension;
                    for (int k = 0; k < mesh_type::dimension; k++) {
                        mm(row + k, col + k) = scalar_matrix(i, j);
                    }
                }
            }

            return mm;
        };

        int offset = 0;
        for (size_t i = 0; i < fcs.size(); i++) {
            const auto fdi = faces_infos[i];

            if (fdi.hasUnknowns()) {
                const auto fc = fcs[i];
                const auto facdeg = fdi.degree();
                const auto hF = diameter(msh, fc);
                const auto fb = make_scalar_monomial_basis(msh, fc, facdeg);
                const auto fbs =
                    vector_basis_size(facdeg, mesh_type::dimension - 1, mesh_type::dimension);

                const matrix_type mass_F = make_mass_matrix(msh, fc, fb);

                mm.block(offset, offset, fbs, fbs) += (1.0 / hF) * to_vector(mass_F);

                offset += fbs;
            }
        }
        assert(offset == num_faces_dofs);

        return mm;
    }

  public:
    SplittingIteration(const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                       const MeshDegreeInfo<mesh_type> &degree_infos,
                       const TimeStep<scalar_type> &current_step)
        : m_verbose(rp.m_verbose), m_time_step(current_step), m_dyna(rp) {

        m_assembler = assembler_type(msh, degree_infos, bnd);
    }

    bool verbose(void) const { return m_verbose; }

    void verbose(bool v) { m_verbose = v; }

    void initialize(const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                    const MeshDegreeInfo<mesh_type> &degree_infos,
                    const std::vector<matrix_type> &stab_precomputed,
                    const StabCoeffManager<scalar_type> &stab_manager,
                    MultiTimeField<scalar_type> &fields) {
        m_dyna.prediction(msh, degree_infos, m_time_step, fields);

        vector_mechanics_hho_assembler assembler(msh, degree_infos, bnd);

        const bool mixed_order = rp.m_cell_degree > rp.m_face_degree;

        if (!mixed_order) {
            m_stab_zFF.reserve(msh.cells_size());
        }

        m_lhs.reserve(msh.cells_size());
        m_stab_uT.reserve(msh.cells_size());
        auto depl = fields.getCurrentField(FieldName::DEPL_CELLS);

        for (auto &cl : msh) {
            const auto cell_i = msh.lookup(cl);

            const auto cell_infos = degree_infos.cellDegreeInfo(msh, cl);
            const auto num_cell_dofs = vector_cell_dofs(msh, cell_infos);

            const auto faces_infos = cell_infos.facesDegreeInfo();
            const auto num_faces_dofs = vector_faces_dofs(msh, faces_infos);

            const auto beta_s = stab_manager.getValue(msh, cl);

            matrix_type stab = beta_s * _stab(msh, cl, rp, degree_infos, stab_precomputed);

            const vector_type Sft_uT =
                stab.bottomLeftCorner(num_faces_dofs, num_cell_dofs) * depl.at(cell_i);
            m_stab_uT.push_back(Sft_uT);

            matrix_type lhs;

            if (mixed_order) {
                lhs = stab.bottomRightCorner(num_faces_dofs, num_faces_dofs);
            } else {
                lhs = beta_s * _mass_term(msh, cl, degree_infos);
                const matrix_type zFF =
                    stab.bottomRightCorner(num_faces_dofs, num_faces_dofs) - lhs;
                m_stab_zFF.push_back(zFF);
            }
            m_lhs.push_back(lhs);

            const vector_type rhs = vector_type::Zero(num_faces_dofs);

            assembler.assemble(msh, cl, bnd, lhs, rhs);
        }

        assembler.finalize();

        m_solver.analyzePattern(assembler.LHS);
        m_solver.factorize(assembler.LHS);

        m_system_displ = assembler.RHS;

        if (m_solver.info() != Eigen::Success) {
            throw std::runtime_error("Fail to solve the linear solver.");
        }
    }

    template <typename LoadFunction>
    AssemblyInfo assemble(const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                          const MeshDegreeInfo<mesh_type> &degree_infos, const LoadFunction &lf,
                          const std::vector<matrix_type> &gradient_precomputed,
                          const std::vector<matrix_type> &stab_precomputed, behavior_type &behavior,
                          StabCoeffManager<scalar_type> &stab_manager,
                          MultiTimeField<scalar_type> &fields) {
        elem_type elem;
        AssemblyInfo ai;

        // set RHS to zero
        m_assembler.initialize();

        const bool small_def = (behavior.getDeformation() == SMALL_DEF);

        const bool mixed_order = rp.m_cell_degree > rp.m_face_degree;

        // Like if it is an implicit scheme
        auto current_time = m_time_step.end_time();
        auto depl = fields.getCurrentField(FieldName::DEPL);

        std::vector<vector_type> resi_cells;
        resi_cells.reserve(msh.cells_size());

        auto rlf = [&lf, &current_time](const point<scalar_type, mesh_type::dimension> &p) -> auto {
            return lf(p, current_time);
        };

        timecounter tc, ttot;

        ttot.tic();

        for (auto &cl : msh) {
            const auto cell_i = msh.lookup(cl);

            const auto huT = depl.at(cell_i);

            const auto cell_infos = degree_infos.cellDegreeInfo(msh, cl);
            const auto num_cell_dofs = vector_cell_dofs(msh, cell_infos);

            const auto faces_infos = cell_infos.facesDegreeInfo();
            const auto num_faces_dofs = vector_faces_dofs(msh, faces_infos);

            // Gradient Reconstruction
            // std::cout << "Grad" << std::endl;
            tc.tic();
            matrix_type GT = _gradrec(msh, cl, rp, degree_infos, small_def, gradient_precomputed);
            tc.toc();
            ai.m_time_gradrec += tc.elapsed();

            // Mechanical Computation

            tc.tic();
            // std::cout << "Elem" << std::endl;
            elem.compute(msh, cl, bnd, rp, degree_infos, rlf, GT, huT, m_time_step, behavior,
                         stab_manager, small_def, false);

            vector_type rhs = elem.RTF.tail(num_faces_dofs);
            resi_cells.push_back(elem.RTF.head(num_cell_dofs));

            rhs -= m_stab_uT.at(cell_i);
            if (!mixed_order) {
                rhs -= m_stab_zFF[cell_i] * huT.tail(num_faces_dofs);
            }

            tc.toc();
            ai.m_time_elem += tc.elapsed();
            ai.m_time_law += elem.time_law;

            m_assembler.assemble_rhs(msh, cl, bnd, m_lhs.at(cell_i), rhs);
        }

        fields.setCurrentField(FieldName::RESI_CELLS, resi_cells);

        m_assembler.impose_neumann_boundary_conditions(msh, bnd);
        m_assembler.finalize();

        ttot.toc();
        ai.m_time_assembly = ttot.elapsed();
        ai.m_linear_system_size = m_assembler.LHS.rows();
        return ai;
    }

    SolveInfo solve(const solvers::LinearSolverType &type) {
        timecounter tc;

        // std::cout << "RHS" << m_assembler.RHS.transpose() << std::endl;

        m_system_displ_prev = m_system_displ;

        tc.tic();
        m_system_displ = m_solver.solve(m_assembler.RHS);
        tc.toc();

        // std::cout << "SOL" << m_system_displ.transpose() << std::endl;

        return SolveInfo(m_assembler.LHS.rows(), m_assembler.LHS.nonZeros(), tc.elapsed());
    }

    scalar_type postprocess(const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                            const MeshDegreeInfo<mesh_type> &degree_infos,
                            MultiTimeField<scalar_type> &fields) {
        timecounter tc;
        tc.tic();

        auto depl = fields.getCurrentField(FieldName::DEPL);

        // Update cell
        for (auto &cl : msh) {
            const auto cell_i = msh.lookup(cl);

            const vector_type udT = m_assembler.take_local_solution(msh, cl, bnd, m_system_displ);

            // Update element U^{i+1} = U^i + delta U^i
            depl.at(cell_i).tail(udT.size()) = udT;

            // std::cout << depl.at(cell_i).transpose() << std::endl;
        }
        fields.setCurrentField(FieldName::DEPL, depl);

        // Update  unknowns
        // Update face Uf^{i+1} = Uf^i + delta Uf^i

        std::vector<vector_type> depl_faces;
        depl_faces.reserve(msh.faces_size());
        for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++) {
            const auto fc = *itor;
            depl_faces.push_back(m_assembler.take_local_solution(msh, fc, bnd, m_system_displ));
        }
        fields.setCurrentField(FieldName::DEPL_FACES, depl_faces);

        tc.toc();
        return tc.elapsed();
    }

    bool convergence(const param_type &rp, const size_t iter,
                     const MultiTimeField<scalar_type> &fields) {

        auto depl_faces = fields.getCurrentField(FieldName::DEPL_FACES);

        // norm of the solution
        auto error_un = norm(depl_faces);

        if (error_un <= scalar_type(10E-15)) {
            error_un = scalar_type(10E16);
        }

        // norm of the rhs
        const vector_type incr = m_system_displ - m_system_displ_prev;
        scalar_type max_error = 0.0;
        for (size_t i = 0; i < incr.size(); i++)
            max_error = std::max(max_error, std::abs(incr(i)));

        // norm of the increment
        const scalar_type error_incr = incr.norm();
        scalar_type relative_displ = error_incr / error_un;

        if (iter == 0) {
            relative_displ = 1;
        }

        if (m_verbose) {
            std::string s_iter = "   " + std::to_string(iter) + "               ";
            s_iter.resize(9);

            if (iter == 0) {
                std::cout << "--------------------------------------------------------------"
                          << std::endl;
                std::cout << "| Iteration | Norme l2 incr | Relative incr  | Maximum error |"
                          << std::endl;
                std::cout << "--------------------------------------------------------------"
                          << std::endl;
            }
            std::ios::fmtflags f(std::cout.flags());
            std::cout.precision(5);
            std::cout.setf(std::iostream::scientific, std::iostream::floatfield);
            std::cout << "| " << s_iter << " |   " << error_incr << " |   " << relative_displ
                      << "  |  " << max_error << "  |" << std::endl;
            std::cout << "--------------------------------------------------------------"
                      << std::endl;
            std::cout.flags(f);
        }

        const scalar_type error = relative_displ;

        if (!isfinite(error))
            throw std::runtime_error("Norm of residual is not finite");

        if (error <= rp.getConvergenceCriteria()) {
            return true;
        } else {
            return false;
        }
    }

    scalar_type post_convergence(const mesh_type &msh, const bnd_type &bnd, const param_type &rp,
                                 const MeshDegreeInfo<mesh_type> &degree_infos,
                                 const std::vector<matrix_type> &stab_precomputed,
                                 const StabCoeffManager<scalar_type> &stab_manager,
                                 MultiTimeField<scalar_type> &fields) {
        timecounter tc;
        tc.tic();

        const auto depl = fields.getCurrentField(FieldName::DEPL);
        auto resi_cells = fields.getCurrentField(FieldName::RESI_CELLS);

        // Update cell
        for (auto &cl : msh) {
            const auto cell_i = msh.lookup(cl);

            const auto cell_infos = degree_infos.cellDegreeInfo(msh, cl);
            const auto num_cell_dofs = vector_cell_dofs(msh, cell_infos);

            const auto beta_s = stab_manager.getValue(msh, cl);
            const matrix_type stab = beta_s * _stab(msh, cl, rp, degree_infos, stab_precomputed);

            // Update residual
            resi_cells[cell_i] -= stab.topLeftCorner(num_cell_dofs, stab.cols()) * depl[cell_i];

            // std::cout << depl.at(cell_i).transpose() << std::endl;
        }
        fields.setCurrentField(FieldName::RESI_CELLS, resi_cells);

        tc.toc();
        return tc.elapsed();
    }
};
} // namespace mechanics
} // namespace disk