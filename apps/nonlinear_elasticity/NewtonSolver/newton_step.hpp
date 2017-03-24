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

 // NewtonRaphson_step

#pragma once

#include <iostream>

#include <sstream>
#include <string>


#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "hho/hho.hpp"
#include "../hho_nl.hpp"

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
class NewtonRaphson_step_elasticity
{
    typedef Mesh                                       mesh_type;
    typedef typename mesh_type::scalar_type            scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;

    typedef disk::scaled_monomial_vector_basis<mesh_type, cell_type>    cell_basis_type;
    typedef disk::scaled_monomial_vector_basis<mesh_type, face_type>    face_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    div_cell_basis_type;

    typedef dynamic_matrix<scalar_type>         matrix_dynamic;
    typedef dynamic_vector<scalar_type>         vector_dynamic;


    typedef disk::gradient_reconstruction_elas<  mesh_type, cell_basis_type, cell_quadrature_type,
                                                face_basis_type, face_quadrature_type>                          grad_type;

    typedef disk::sgradient_reconstruction_elas<  mesh_type, cell_basis_type, cell_quadrature_type,
                                                face_basis_type, face_quadrature_type>                          sgrad_type;

    typedef disk::elas_like_stabilization<  mesh_type, cell_basis_type, cell_quadrature_type,
                                            face_basis_type, face_quadrature_type>                              stab_type;
    typedef disk::diffusion_like_static_condensation<   mesh_type, cell_basis_type, cell_quadrature_type,
                                                        face_basis_type, face_quadrature_type>                  statcond_type;
    typedef disk::assembler_nl_elas<   mesh_type, face_basis_type, face_quadrature_type>                        assembler_type;

    typedef disk::divergence_reconstruction_elas<mesh_type, cell_basis_type, cell_quadrature_type,
                                          face_basis_type, face_quadrature_type, div_cell_basis_type,
                                          cell_quadrature_type>                                                 div_type;

    size_t m_cell_degree, m_face_degree, m_degree;

    typedef typename assembler_type::sparse_matrix_type       sparse_matrix_type;
    typedef typename assembler_type::vector_type            vector_type;

    #ifdef HAVE_INTEL_MKL
        typedef Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver_type;
    #else
        typedef Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver_type;
    #endif

    sparse_matrix_type     m_system_matrix;
    vector_type            m_system_rhs, m_system_solution;

    solver_type             solver;

    const mesh_type& m_msh;

    std::vector<vector_dynamic>        m_postprocess_data, m_solution_data;
    std::vector<vector_dynamic>         m_solution_cells, m_solution_faces, m_solution_lagr;

    bool m_verbose;



solver_info
    solve_computelu(const bool reactualisze_next)
    {
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


        //free_memomry (pas très beau)
        m_system_matrix = sparse_matrix_type(0,0);
        m_system_rhs = vector_type::Zero(0);
        if (reactualisze_next) {
           sparse_matrix_type id = sparse_matrix_type(1,1);
           id.setIdentity();
           solver.analyzePattern(id);
           solver.factorize(id);
        }

        return si;
    }

solver_info
    solve_reuselu(const bool reactualisze_next)
    {
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
         m_system_solution = solver.solve(m_system_rhs);
        tc.toc();
        si.time_solver = tc.to_double();


        //free_memomry
        m_system_matrix = sparse_matrix_type(0,0);
        m_system_rhs = vector_type::Zero(0);
        if (reactualisze_next) {
           sparse_matrix_type id = sparse_matrix_type(1,1);
           id.setIdentity();
           solver.analyzePattern(id);
           solver.factorize(id);
        }

        return si;
    }

public:
    NewtonRaphson_step_elasticity(const mesh_type& msh, size_t degree, int l = 0)
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
        m_degree = degree;

    }

    bool    verbose(void) const     { return m_verbose; }
    void    verbose(bool v)         { m_verbose = v; }



    void
    initialize( const std::vector<vector_dynamic>& initial_solution_cells,
                const std::vector<vector_dynamic>& initial_solution_faces,
                const std::vector<vector_dynamic>& initial_solution_lagr,
                const std::vector<vector_dynamic>& initial_solution)
    {
        assert(m_msh.cells_size() == initial_solution_cells.size());
        m_solution_cells = initial_solution_cells;
        assert(m_msh.cells_size() == m_solution_cells.size());

        assert(m_msh.faces_size() == initial_solution_faces.size());
        m_solution_faces = initial_solution_faces;
        assert(m_msh.faces_size() == m_solution_faces.size());

        assert(m_msh.boundary_faces_size() == initial_solution_lagr.size());
        m_solution_lagr = initial_solution_lagr;
        assert(m_msh.boundary_faces_size() == m_solution_lagr.size());

        assert(m_msh.cells_size() == initial_solution.size());
        m_solution_data = initial_solution;
        assert(m_msh.cells_size() == m_solution_data.size());

    }

    template<typename LoadFunction, typename BoundaryConditionFunction>
    assembly_info
    assemblelin(const LoadFunction& lf, const BoundaryConditionFunction& bcf)
    {


       sgrad_type gradrec(m_degree);
       stab_type stab(m_degree);
       //div_type   divrec(m_degree);

       statcond_type statcond(m_degree);

       assembler_type assembler(m_msh, m_degree);


       assembly_info ai;
       bzero(&ai, sizeof(ai));

       timecounter tc;

       scalar_type mu      = 1.0;
       scalar_type lambda  = 0.0;
       const size_t DIM = m_msh.dimension;

       size_t i = 0;

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

          //divrec.compute(m_msh, cl);
          tc.tic();

          std::cout << "SOLUTIONPRED" << std::endl;
          std::cout << m_solution_data.at(i) << std::endl;

          std::cout << "gradrec" << std::endl;
          std::cout << gradrec.oper << std::endl;
          std::cout << "stab" << std::endl;
          std::cout << stab.data << std::endl;

          /////// LINEAIRE ////////////////


                  dynamic_matrix<scalar_type> lhs = 2.0*mu*gradrec.data + 2.0*mu*stab.data;// + lambda * divrec.data;


                  dynamic_vector<scalar_type> rhsext = disk::compute_rhs_ext<cell_basis_type, cell_quadrature_type>
                                                      (m_msh, cl, lf, m_cell_degree,  m_solution_data.at(i));

                  dynamic_vector<scalar_type> rhs = -( lhs * m_solution_data.at(i) - rhsext);

                                                      std::cout << "rhs" << std::endl;
                                                      std::cout << rhs << std::endl;

          auto scnp = statcond.compute(m_msh, cl, lhs, rhs, true);
          tc.toc();
          ai.time_statcond += tc.to_double();

          std::cout << "mat_cond" << std::endl;
          std::cout << scnp.first << std::endl;
          std::cout << "rhs_cond" << std::endl;
          std::cout << scnp.second << std::endl;

          assembler.assemble(m_msh, cl, scnp);
          i++;
       }

       assembler.impose_boundary_conditions(m_msh, bcf, m_solution_faces, m_solution_lagr);
       assembler.finalize(m_system_matrix, m_system_rhs);

       ai.linear_system_size = m_system_matrix.rows();
       return ai;
    }

    template<typename LoadFunction, typename BoundaryConditionFunction>
    assembly_info
    assemble(const LoadFunction& lf, const BoundaryConditionFunction& bcf)
    {


        grad_type gradrec(m_degree);
        stab_type stab(m_degree);

        statcond_type statcond(m_degree);

        assembler_type assembler(m_msh, m_degree);


        assembly_info ai;
        bzero(&ai, sizeof(ai));

        timecounter tc;

        scalar_type mu      = 1.0;
        scalar_type lambda  = 0.0;
        const size_t DIM = m_msh.dimension;

        size_t i = 0;

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


            /////// NON LINEAIRE /////////
            auto gtu = gradrec.oper * m_solution_data.at(i);
            auto stu = stab.data * m_solution_data.at(i);

            auto elem = disk::compute_elem<cell_basis_type, cell_quadrature_type, mesh_type>(m_msh, cl, gtu, m_degree);

            dynamic_vector<scalar_type> rhs_ext = disk::compute_rhs_ext<cell_basis_type, cell_quadrature_type>
                                                (m_msh, cl, lf, m_cell_degree, m_solution_data.at(i));

            dynamic_matrix<scalar_type> lhs = disk::assemble_lhs<scalar_type>(gradrec.oper, stab.data, elem.first, 2.0*mu);

            dynamic_vector<scalar_type> rhs = - disk::assemble_rhs<scalar_type>(gradrec.oper, stu, elem.second, rhs_ext, 2.0*mu);

            /////// LINEAIRE ////////////////

/*
            dynamic_matrix<scalar_type> lhs = 2.0*mu*gradrec.data + 2.0*mu*stab.data + lambda * divrec.data;


            dynamic_vector<scalar_type> rhs = disk::compute_rhs_ext<cell_basis_type, cell_quadrature_type>
                                                (m_msh, cl, lf, m_cell_degree,  m_solution_data.at(i));*/

            auto scnp = statcond.compute(m_msh, cl, lhs, rhs, true);
            tc.toc();
            ai.time_statcond += tc.to_double();

            assembler.assemble(m_msh, cl, scnp);



            std::cout << "SOLUTIONPRED" << std::endl;
            std::cout << m_solution_data.at(i) << std::endl;

            std::cout << "Gu" << std::endl;
            std::cout << gtu << std::endl;
            std::cout << "STu" << std::endl;
            std::cout << stu << std::endl;


            std::cout << "mat_cond" << std::endl;
            std::cout << scnp.first << std::endl;
            std::cout << "rhs_cond" << std::endl;
            std::cout << scnp.second << std::endl;

            i++;

        }

         assembler.impose_boundary_conditions(m_msh, bcf, m_solution_faces, m_solution_lagr);
         assembler.finalize(m_system_matrix, m_system_rhs);

         ai.linear_system_size = m_system_matrix.rows();
        return ai;
    }


//         template<typename LoadFunction, typename BoundaryConditionFunction>
//     assembly_info
//     assemble(const LoadFunction& lf, const BoundaryConditionFunction& bcf, const std::vector<matrix_dynamic>& offline_data)
//     {
//         assert(offline_data.size() == m_msh.cells_size());
//         statcond_type statcond(m_degree);
//
//         assembler_type assembler(m_msh, m_degree);
//
//         assembly_info ai;
//         bzero(&ai, sizeof(ai));
//
//         timecounter tc;
//
//         scalar_type mu      = 1.0;
//         scalar_type lambda  = 1.0;
//         const size_t DIM = m_msh.dimension;
//
//         size_t i = 0;
//         for (auto& cl : m_msh)
//         {
//             tc.tic();
//
//             dynamic_vector<scalar_type> rhs = disk::compute_rhs_diffusion<cell_basis_type, cell_quadrature_type>
//                                                 (m_msh, cl, lf, m_cell_degree, offline_data.at(i),  m_solution_data.at(i));
//             auto scnp = statcond.compute(m_msh, cl, offline_data.at(i), -rhs, true);
//             tc.toc();
//             ai.time_statcond += tc.to_double();
//
//             assembler.assemble(m_msh, cl, scnp);
//
//             i++;
//         }
//
//          assembler.impose_boundary_conditions(m_msh, bcf, m_solution_faces, m_solution_lagr);
//          assembler.finalize(m_system_matrix, m_system_rhs);
//
//          ai.linear_system_size = m_system_matrix.rows();
//         return ai;
//     }

    solver_info
    solve(const bool reactualize,  const bool reactualisze_next)
    {
        if (reactualize)
            return solve_computelu(reactualisze_next);
        else
            return solve_reuselu(reactualisze_next);
    }

    template<typename LoadFunction>
    postprocess_info
    postprocesslin(const LoadFunction& lf)
    {
       sgrad_type gradrec(m_degree);
       stab_type stab(m_degree);
       //div_type   divrec(m_degree);

       statcond_type statcond(m_degree);

       face_basis_type face_basis(m_degree);
       size_t fbs = face_basis.size();
       cell_basis_type cell_basis(m_degree);
       size_t cbs = cell_basis.size();

       postprocess_info pi;

       m_postprocess_data.clear();

       m_postprocess_data.reserve(m_msh.cells_size());

       timecounter tc;
       tc.tic();

       scalar_type lambda = 1.0;
       scalar_type mu = 1.0;

       size_t i =0;
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
          //divrec.compute(m_msh, cl);


          /////// LINEAIRE ////////////////


                  dynamic_matrix<scalar_type> lhs = 2.0*mu*gradrec.data + 2.0*mu*stab.data;// + lambda * divrec.data;


                  dynamic_vector<scalar_type> rhsext = disk::compute_rhs_ext<cell_basis_type, cell_quadrature_type>
                                                      (m_msh, cl, lf, m_cell_degree,  m_solution_data.at(i));

          dynamic_vector<scalar_type> rhs = -( lhs * m_solution_data.at(i) - rhsext);
          dynamic_vector<scalar_type> rhs_cell = rhs.block(0,0, cbs, 1);
          dynamic_vector<scalar_type> x = statcond.recover(m_msh, cl, lhs, rhs_cell, xFs);
          m_postprocess_data.push_back(x);

          std::cout << "sol" << std::endl;
          std::cout << xFs << std::endl;
          std::cout << "recover" << std::endl;
          std::cout << x << std::endl;
       }
       tc.toc();

       pi.time_postprocess = tc.to_double();

       return pi;
    }

    template<typename LoadFunction>
    postprocess_info
    postprocess(const LoadFunction& lf)
    {
        grad_type gradrec(m_degree);
        stab_type stab(m_degree);

        statcond_type statcond(m_degree);

        face_basis_type face_basis(m_degree);
        size_t fbs = face_basis.size();
        cell_basis_type cell_basis(m_degree);
        size_t cbs = cell_basis.size();

        postprocess_info pi;

        m_postprocess_data.clear();

        m_postprocess_data.reserve(m_msh.cells_size());

        timecounter tc;
        tc.tic();

        scalar_type lambda = 1.0;
        scalar_type mu =1.0;

        size_t i =0;
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

            /////// NON LINEAIRE /////////
            auto gtu = gradrec.oper * m_solution_data.at(i);
            auto stu = stab.data * m_solution_data.at(i);

            auto elem = disk::compute_elem<cell_basis_type, cell_quadrature_type, mesh_type>(m_msh, cl, gtu, m_degree);

            dynamic_vector<scalar_type> rhs_ext = disk::compute_rhs_ext<cell_basis_type, cell_quadrature_type>
                                                (m_msh, cl, lf, m_cell_degree, m_solution_data.at(i));

            dynamic_matrix<scalar_type> lhs = disk::assemble_lhs<scalar_type>(gradrec.oper, stab.data, elem.first, 2.0*mu);

            dynamic_vector<scalar_type> rhs = - disk::assemble_rhs<scalar_type>(gradrec.oper, stu, elem.second, rhs_ext, 2.0*mu);

            /////// LINEAIRE ////////////////

/*
            dynamic_matrix<scalar_type> lhs = 2.0*mu*gradrec.data + 2.0*mu*stab.data + lambda * divrec.data;


            dynamic_vector<scalar_type> rhs = disk::compute_rhs_ext<cell_basis_type, cell_quadrature_type>
                                                (m_msh, cl, lf, m_cell_degree,  m_solution_data.at(i));*/

            dynamic_vector<scalar_type> rhs_cell = rhs.block(0,0, cbs, 1);
            dynamic_vector<scalar_type> x = statcond.recover(m_msh, cl, lhs, rhs_cell, xFs);
            m_postprocess_data.push_back(x);
        }
        tc.toc();

        pi.time_postprocess = tc.to_double();

        return pi;
    }


//     template<typename LoadFunction>
//     postprocess_info
//     postprocess(const LoadFunction& lf,  const std::vector<matrix_dynamic>& offline_data)
//     {
//         assert(offline_data.size() == m_msh.cells_size());
//         statcond_type statcond(m_degree);
//
//         face_basis_type face_basis(m_degree);
//         size_t fbs = face_basis.size();
//         cell_basis_type cell_basis(m_degree);
//         size_t cbs = cell_basis.size();
//
//         postprocess_info pi;
//
//         m_postprocess_data.clear();
//
//         m_postprocess_data.reserve(m_msh.cells_size());
//
//         timecounter tc;
//         tc.tic();
//
//         size_t i =0;
//         for (auto& cl : m_msh)
//         {
//             auto fcs = faces(m_msh, cl);
//             auto num_faces = fcs.size();
//
//             dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(num_faces*fbs);
//
//             for (size_t face_i = 0; face_i < num_faces; face_i++)
//             {
//                 auto fc = fcs[face_i];
//                 auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), fc);
//                 if (!eid.first)
//                     throw std::invalid_argument("This is a bug: face not found");
//
//                 auto face_id = eid.second;
//
//                 dynamic_vector<scalar_type> xF = dynamic_vector<scalar_type>::Zero(fbs);
//                 xF = m_system_solution.block(face_id * fbs, 0, fbs, 1);
//                 xFs.block(face_i * fbs, 0, fbs, 1) = xF;
//             }
//
//             dynamic_vector<scalar_type> rhs = disk::compute_rhs_diffusion<cell_basis_type, cell_quadrature_type>
//                                                 (m_msh, cl, lf, m_cell_degree, offline_data.at(i),  m_solution_data.at(i));
//             dynamic_vector<scalar_type> rhs_cell = -rhs.block(0,0, cbs, 1);
//             dynamic_vector<scalar_type> x = statcond.recover(m_msh, cl,  offline_data.at(i), rhs_cell, xFs);
//             m_postprocess_data.push_back(x);
//
//             i++;
//         }
//         tc.toc();
//
//         pi.time_postprocess = tc.to_double();
//
//         return pi;
//     }

    void
    update_solution()
    {
        assert(m_postprocess_data.size() == m_solution_data.size());

        for(size_t i=0; i < m_postprocess_data.size(); i++){
            assert(m_solution_data.at(i).size() == m_postprocess_data.at(i).size());
            m_solution_data.at(i) += m_postprocess_data.at(i);
        }

        cell_basis_type cell_basis(m_degree);
        size_t cbs = cell_basis.size();

        assert(m_postprocess_data.size() == m_solution_cells.size());

        for(size_t i=0; i < m_postprocess_data.size(); i++){
            assert(m_solution_cells.at(i).size() == cbs);
            m_solution_cells.at(i) += (m_postprocess_data.at(i)).block(0,0,cbs,1);
        }

        face_basis_type face_basis(m_degree);
        size_t fbs = face_basis.size();


        for(size_t i=0; i < m_solution_faces.size(); i++){
            assert(m_solution_faces.at(i).size() == fbs);
            m_solution_faces.at(i) += m_system_solution.block(i * fbs, 0, fbs, 1);
        }

        const size_t lagrange_offset = m_solution_faces.size() * fbs;
        for(size_t i=0; i < m_solution_lagr.size(); i++){
            assert(m_solution_lagr.at(i).size() == fbs);
            m_solution_lagr.at(i) += m_system_solution.block(lagrange_offset + i * fbs, 0, fbs, 1);
        }
    }

    bool test_convergence(const scalar_type epsilon, const size_t iter)
    {
      // a calculer erreur
      scalar_type relative_error(0);
      scalar_type max_error(0);

      relative_error = m_system_rhs.dot(m_system_rhs);

      std::cout << "m_rhs_systeme" << std::endl;
      std::cout << m_system_rhs << std::endl;

      for (size_t j = 0; j < m_system_rhs.size(); j++){
         scalar_type test_error = std::abs(m_system_rhs(j));
         if( test_error > max_error)
              max_error = test_error;
      }

// OLD TEST based on increment
//         for(size_t i=0; i < m_postprocess_data.size(); i++ ){
//             relative_error += (m_postprocess_data.at(i)).dot(m_postprocess_data.at(i));
//             relative_error = m_system_rhs.dot(m_system_rhs);
//
//             for (size_t j = 0; j < (m_postprocess_data.at(i)).size(); j++){
//                 scalar_type test_error = std::abs((m_postprocess_data.at(i))(j));
//                 if( test_error > max_error)
//             }
//         }

      relative_error = sqrt(relative_error);

      std::string s_iter = "   " + std::to_string(iter) + "               ";
      s_iter.resize(9);

      if(iter == 0){
         std::cout << "----------------------------------------------" << std::endl;
         std::cout << "| Iteration | Relative error | Maximum error |" << std::endl;
         std::cout << "----------------------------------------------" << std::endl;

      }
      std::ios::fmtflags f( std::cout.flags() );
      std::cout.precision(5);
      std::cout.setf(std::iostream::scientific, std::iostream::floatfield);
      std::cout << "| " << s_iter << " |   " << relative_error << "  |  " << max_error << "  |" << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
      std::cout.flags( f );

      if(relative_error <= epsilon){
         return true;
      }
      else {
         return false;
      }

   }

   void
    save_solutions( std::vector<vector_dynamic>& solution_cells,
                    std::vector<vector_dynamic>& solution_faces,
                    std::vector<vector_dynamic>& solution_lagr,
                    std::vector<vector_dynamic>& solution)
    {
        assert(m_solution_cells.size() == solution_cells.size());
        solution_cells = m_solution_cells;

        assert(m_solution_faces.size() == solution_faces.size());
        solution_faces = m_solution_faces;

        assert(m_solution_lagr.size() == solution_lagr.size());
        solution_lagr = m_solution_lagr;

        assert(m_solution_data.size() == solution.size());
        solution = m_solution_data;
    }
};
