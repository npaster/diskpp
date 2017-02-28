/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016, 2017 - matteo.cicuttin@enpc.fr
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

#include <sstream>

#include "../../config.h"

#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "hho/hho.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

struct assembly_info
{
    size_t  linear_system_size;
    double  time_gradrec, time_statcond, time_stab;
};

struct solver_info
{
    double  time_solver;
};

struct postprocess_info
{
    double  time_postprocess;
};

template<typename Mesh>
class diffusion_solver
{
    typedef Mesh                                       mesh_type;
    typedef typename mesh_type::scalar_type            scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;

    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

    typedef
    disk::basis_quadrature_data<mesh_type,
                                disk::scaled_monomial_scalar_basis,
                                disk::quadrature> bqdata_type;

    typedef disk::gradient_reconstruction_bq<bqdata_type>               gradrec_type;
    typedef disk::diffusion_like_stabilization_bq<bqdata_type>          stab_type;
    typedef disk::diffusion_like_static_condensation_bq<bqdata_type>    statcond_type;
    typedef disk::assembler<mesh_type, face_basis_type, face_quadrature_type> assembler_type;

    size_t m_cell_degree, m_face_degree;

    bqdata_type     m_bqd;

    typename assembler_type::sparse_matrix_type     m_system_matrix;
    typename assembler_type::vector_type            m_system_rhs, m_system_solution;

    const mesh_type& m_msh;

    std::vector<dynamic_vector<scalar_type>>        m_postprocess_data;

    bool m_verbose;

public:
    diffusion_solver(const mesh_type& msh, size_t degree, int l = 0)
        : m_msh(msh), m_verbose(false)
    {
        if ( l < -1 or l > 1)
        {
            std::cout << "'l' should be -1, 0 or 1. Reverting to 0." << std::endl;
            l = 0;
        }

        if (degree == 0 && l == -1)
        {
            std::cout << "'l' should be 0 or 1. Reverting to 0." << std::endl;
            l = 0;
        }

        m_cell_degree = degree + l;
        m_face_degree = degree;

        m_bqd = bqdata_type(m_cell_degree, m_face_degree);
    }

    bool    verbose(void) const     { return m_verbose; }
    void    verbose(bool v)         { m_verbose = v; }

    template<typename LoadFunction, typename BoundaryConditionFunction>
    assembly_info
    assemble(const LoadFunction& lf, const BoundaryConditionFunction& bcf)
    {
        auto gradrec    = gradrec_type(m_bqd);
        auto stab       = stab_type(m_bqd);
        auto statcond   = statcond_type(m_bqd);
        auto assembler  = assembler_type(m_msh, m_face_degree);

        assembly_info ai;
        bzero(&ai, sizeof(ai));

        timecounter tc;

        for (auto& cl : m_msh)
        {
            tc.tic();
            gradrec.compute(m_msh, cl);
            tc.toc();
            ai.time_gradrec += tc.to_double();

            tc.tic();
            stab.compute(m_msh, cl, gradrec.oper);
            tc.toc();
            ai.time_stab += tc.to_double();

            tc.tic();
            auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(m_msh, cl, lf, m_cell_degree);
            dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
            auto scnp = statcond.compute(m_msh, cl, loc, cell_rhs);
            tc.toc();
            ai.time_statcond += tc.to_double();

            assembler.assemble(m_msh, cl, scnp);
        }

        assembler.impose_boundary_conditions(m_msh, bcf);
        assembler.finalize(m_system_matrix, m_system_rhs);

        ai.linear_system_size = m_system_matrix.rows();
        return ai;
    }

    solver_info
    solve(void)
    {
#ifdef HAVE_INTEL_MKL
        Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
#else
        Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
#endif

        solver_info si;

        size_t systsz = m_system_matrix.rows();
        size_t nnz = m_system_matrix.nonZeros();

        if (verbose())
        {
            std::cout << "Starting linear solver..." << std::endl;
            std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
            std::cout << " * Matrix fill: " << 100.0*double(nnz)/(systsz*systsz) << "%" << std::endl;
        }

        timecounter tc;

        tc.tic();
        solver.analyzePattern(m_system_matrix);
        solver.factorize(m_system_matrix);
        m_system_solution = solver.solve(m_system_rhs);
        tc.toc();
        si.time_solver = tc.to_double();

        return si;
    }

    template<typename LoadFunction>
    postprocess_info
    postprocess(const LoadFunction& lf)
    {
        auto gradrec    = gradrec_type(m_bqd);
        auto stab       = stab_type(m_bqd);
        auto statcond   = statcond_type(m_bqd);

        size_t fbs = m_bqd.face_basis.size();

        postprocess_info pi;

        m_postprocess_data.reserve(m_msh.cells_size());

        timecounter tc;
        tc.tic();
        for (auto& cl : m_msh)
        {
            auto fcs = faces(m_msh, cl);
            auto num_faces = fcs.size();

            dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(num_faces*fbs);

            for (size_t face_i = 0; face_i < num_faces; face_i++)
            {
                auto fc = fcs[face_i];
                auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), fc);
                if (!eid.first)
                    throw std::invalid_argument("This is a bug: face not found");

                auto face_id = eid.second;

                dynamic_vector<scalar_type> xF = dynamic_vector<scalar_type>::Zero(fbs);
                xF = m_system_solution.block(face_id * fbs, 0, fbs, 1);
                xFs.block(face_i * fbs, 0, fbs, 1) = xF;
            }

            gradrec.compute(m_msh, cl);
            stab.compute(m_msh, cl, gradrec.oper);
            dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
            auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(m_msh, cl, lf, m_cell_degree);
            dynamic_vector<scalar_type> x = statcond.recover(m_msh, cl, loc, cell_rhs, xFs);
            m_postprocess_data.push_back(x);
        }
        tc.toc();

        pi.time_postprocess = tc.to_double();

        return pi;
    }

    template<typename AnalyticalSolution>
    scalar_type
    compute_l2_error(const AnalyticalSolution& as)
    {
        scalar_type err_dof = 0.0;

        disk::projector_bq<bqdata_type> projk(m_bqd);

        size_t i = 0;
        for (auto& cl : m_msh)
        {
            auto x = m_postprocess_data.at(i++);
            dynamic_vector<scalar_type> true_dof = projk.compute_cell(m_msh, cl, as);
            dynamic_vector<scalar_type> comp_dof = x.block(0,0,true_dof.size(), 1);
            dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
            err_dof += diff_dof.dot(projk.cell_mm * diff_dof);
        }

        return sqrt(err_dof);
    }

    void
    plot_solution(const std::string& filename)
    {
        std::ofstream ofs(filename);

        size_t cell_i = 0;
        for (auto& cl : m_msh)
        {
            auto x = m_postprocess_data.at(cell_i++);
            auto qps = m_bqd.cell_quadrature.integrate(m_msh, cl);
            for (auto& qp : qps)
            {
                auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, qp.point());

                scalar_type pot = 0.0;
                for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_cell_degree).size(); i++)
                    pot += phi[i] * x(i);

                auto tp = qp.point();
                for (size_t i = 0; i < mesh_type::dimension; i++)
                    ofs << tp[i] << " ";
                ofs << pot << std::endl;
            }
        }

        ofs.close();
    }
};