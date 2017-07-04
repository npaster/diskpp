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
#include <fstream>


#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "hho/hho.hpp"
#include "hho/hho_vector.hpp"
#include "hho/hho_nl.hpp"
#include "hho/hho_nl_vector.hpp"
#include "../ElasticityParameters.hpp"
#include "../hyperelasticity_elementary_computation.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

struct assembly_info
{
    size_t  linear_system_size;
    double  time_gradrec, time_potrec, time_statcond, time_stab, time_elem, time_law, time_adapt_stab;
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
class NewtonRaphson_step_hyperelasticity
{
    typedef Mesh                                       mesh_type;
    typedef typename mesh_type::scalar_type            scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;
    typedef disk::quadrature<mesh_type, cell_type>      grad_quadrature_type;

    typedef disk::scaled_monomial_vector_basis<mesh_type, cell_type>    cell_basis_type;
    typedef disk::scaled_monomial_vector_basis<mesh_type, face_type>    face_basis_type;
    typedef disk::scaled_monomial_matrix_basis<mesh_type, cell_type>    grad_basis_type;

    typedef dynamic_matrix<scalar_type>         matrix_dynamic;
    typedef dynamic_vector<scalar_type>         vector_dynamic;


    typedef disk::gradient_reconstruction_elas_full<  mesh_type, cell_basis_type, cell_quadrature_type,
                                             face_basis_type, face_quadrature_type,
                                           grad_basis_type, grad_quadrature_type>                               gradrec_type;

    typedef disk::gradient_reconstruction_elas<  mesh_type, cell_basis_type, cell_quadrature_type,
                                                            face_basis_type, face_quadrature_type>              deplrec_type;

    typedef Hyperelasticity::Hyperelasticity<  mesh_type, cell_basis_type, cell_quadrature_type,
                                           face_basis_type, face_quadrature_type,
                                           grad_basis_type, grad_quadrature_type>                               hyperelasticity_type;

    typedef disk::elas_like_stabilization_l2<  mesh_type, cell_basis_type, cell_quadrature_type,
                                            face_basis_type, face_quadrature_type>                              stab_type;

   typedef disk::elas_like_stabilization<  mesh_type, cell_basis_type, cell_quadrature_type,
                                                      face_basis_type, face_quadrature_type>                    stab2_type;
    typedef disk::diffusion_like_static_condensation<   mesh_type, cell_basis_type, cell_quadrature_type,
                                                        face_basis_type, face_quadrature_type>                  statcond_type;
    typedef disk::assembler_nl_vector<   mesh_type, face_basis_type, face_quadrature_type>                   assembler_type;


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
    std::vector<scalar_type>           m_adpat_stab;

    bool m_verbose;

    scalar_type m_beta = 1.0;

    ElasticityParameters m_elas_param;

    scalar_type initial_residual;


public:
   NewtonRaphson_step_hyperelasticity(const mesh_type& msh, const size_t degree, const ElasticityParameters elas_param, int l = 0)
   : m_msh(msh), m_verbose(false), m_elas_param(elas_param)
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
        m_beta = elas_param.tau;

    }

    bool    verbose(void) const     { return m_verbose; }
    void    verbose(bool v)         { m_verbose = v; }



    void
    initialize( const std::vector<vector_dynamic>& initial_solution_cells,
                const std::vector<vector_dynamic>& initial_solution_faces,
                const std::vector<vector_dynamic>& initial_solution_lagr,
                const std::vector<vector_dynamic>& initial_solution,
                const std::vector<size_t>& boundary_neumann)
    {
      m_solution_cells.clear();
      m_solution_cells = initial_solution_cells;
      assert(m_msh.cells_size() == m_solution_cells.size());

      m_solution_faces.clear();
      m_solution_faces = initial_solution_faces;
      assert(m_msh.faces_size() == m_solution_faces.size());

      m_solution_lagr.clear();
      m_solution_lagr = initial_solution_lagr;
      assert((m_msh.boundary_faces_size() - number_of_neumann_faces(m_msh, boundary_neumann))
            == m_solution_lagr.size());

      m_solution_data.clear();
      m_solution_data = initial_solution;
      assert(m_msh.cells_size() == m_solution_data.size());
      
      if(m_elas_param.adaptative_stab)
         m_adpat_stab.resize(m_msh.cells_size());
    }

    template<typename LoadFunction, typename BoundaryConditionFunction, typename NeumannFunction>
    assembly_info
    assemble(const LoadFunction& lf, const BoundaryConditionFunction& bcf, const NeumannFunction& g,
            const std::vector<size_t>& boundary_neumann)
    {
        gradrec_type gradrec(m_degree);
        stab_type stab(m_degree);
        hyperelasticity_type hyperelasticity(m_degree);
        deplrec_type deplrec(m_degree);
        //stab2_type stab(m_degree);

        statcond_type statcond(m_degree);

        assembler_type assembler(m_msh, m_degree, boundary_neumann);

        assembly_info ai;
        bzero(&ai, sizeof(ai));

        timecounter tc;

        size_t i = 0;

        for (auto& cl : m_msh)
        {
            tc.tic();
            gradrec.compute(m_msh, cl);
            tc.toc();
            ai.time_gradrec += tc.to_double();

            tc.tic();
            //deplrec.compute(m_msh, cl);
            //stab.compute(m_msh, cl, deplrec.oper);
            stab.compute(m_msh, cl, gradrec.oper);
            tc.toc();
            ai.time_stab += tc.to_double();


            tc.tic();
            hyperelasticity.compute(m_msh, cl, lf, g, boundary_neumann, gradrec.oper, m_solution_data.at(i), m_elas_param, m_elas_param.adaptative_stab);
            //hyperelasticity.compute2(m_msh, cl, lf, g, boundary_neumann, deplrec.oper, m_solution_data.at(i), m_elas_param);
            
            if(m_elas_param.adaptative_stab){
               m_beta = hyperelasticity.beta_adap;
               m_adpat_stab[i] = m_beta;
            }
            else
               m_beta = m_elas_param.tau;

            /////// NON LINEAIRE /////////
            assert( hyperelasticity.K_int.rows() == stab.data.rows());
            assert( hyperelasticity.K_int.cols() == stab.data.cols());
            dynamic_matrix<scalar_type> lhs = hyperelasticity.K_int + m_beta * stab.data;

            assert( hyperelasticity.RTF.rows() == (stab.data * m_solution_data.at(i)).rows());
            assert( hyperelasticity.RTF.cols() == (stab.data * m_solution_data.at(i)).cols());

            dynamic_vector<scalar_type> rhs = hyperelasticity.RTF  - m_beta * stab.data * m_solution_data.at(i) ;

            tc.toc();
            ai.time_elem += tc.to_double();
            ai.time_law += hyperelasticity.time_law;
            ai.time_adapt_stab += hyperelasticity.time_adapt_stab;
            tc.tic();
            auto scnp = statcond.compute(m_msh, cl, lhs, rhs, true);
            tc.toc();
            ai.time_statcond += tc.to_double();

            assembler.assemble(m_msh, cl, scnp);

            i++;
        }

         assembler.impose_boundary_conditions(m_msh, bcf, m_solution_faces, m_solution_lagr, boundary_neumann);
         assembler.finalize(m_system_matrix, m_system_rhs);

         ai.linear_system_size = m_system_matrix.rows();
        return ai;
    }

    void
    saveMatrix(const std::string& filename)
    {
       std::ofstream fichier(filename, std::ios::out | std::ios::trunc);

       for (int k=0; k<m_system_matrix.outerSize(); ++k)
          for (typename Eigen::SparseMatrix<scalar_type>::InnerIterator it(m_system_matrix,k); it; ++it)
         {
            fichier << it.row() << " ; " << it.col() << " ; " << it.value() << std::endl;
         }

       fichier.close();
    }

    void
    test_aurrichio(void)
    {
      face_basis_type face_basis(m_degree);
      const size_t fbs = face_basis.size();
      const size_t nb_faces_dof = fbs * m_msh.faces_size();

      sparse_matrix_type mat_aurri = m_system_matrix.block(0, 0, nb_faces_dof, nb_faces_dof);

      Eigen::SelfAdjointEigenSolver<sparse_matrix_type> es;
      es.compute(mat_aurri);

      if(es.info() != Eigen::Success) {
          std::cerr << "ERROR: Could not compute eigenvalues of the matrix" << std::endl;
      }

      const scalar_type min_ev = es.eigenvalues().minCoeff();
      const scalar_type max_ev = es.eigenvalues().maxCoeff();

      std::cout << "******* Eigenvalues test ********" << std::endl;
      std::cout << "Number of eigenvalues: " << es.eigenvalues().size()  << std::endl;
      std::cout << "Maximum eigenvalue: " << max_ev << std::endl;
      std::cout << "Minimum eigenvalue: " << min_ev << std::endl;
      std::cout << "Conditionning number: " << std::abs(max_ev/min_ev) << std::endl;

      if(min_ev <= -1.0) {
         std::cerr << "The matrix is no mre definite positive" << min_ev << std::endl;
      }
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

//        size_t systsz = m_system_matrix.rows();
//        size_t nnz = m_system_matrix.nonZeros();

       timecounter tc;

       tc.tic();
       solver.analyzePattern(m_system_matrix);
       solver.factorize(m_system_matrix);

       if(solver.info() != Eigen::Success) {
          std::cerr << "ERROR: Could not factorize the matrix" << std::endl;
       }

       m_system_solution = solver.solve(m_system_rhs);
       if(solver.info() != Eigen::Success) {
          std::cerr << "ERROR: Could not solve the linear system" << std::endl;
       }
       tc.toc();
       si.time_solver = tc.to_double();

       return si;
    }



    template<typename LoadFunction, typename NeumannFunction>
    postprocess_info
    postprocess(const LoadFunction& lf, const NeumannFunction& g, const std::vector<size_t>& boundary_neumann)
    {
        gradrec_type gradrec(m_degree);
        stab_type stab(m_degree);
        hyperelasticity_type hyperelasticity(m_degree);
         deplrec_type deplrec(m_degree);
         //stab2_type stab(m_degree);

        statcond_type statcond(m_degree);

        face_basis_type face_basis(m_degree);
        const size_t fbs = face_basis.size();
        cell_basis_type cell_basis(m_degree);
        const size_t cbs = cell_basis.size();

        postprocess_info pi;

        m_postprocess_data.clear();

        m_postprocess_data.reserve(m_msh.cells_size());

        timecounter tc;
        tc.tic();


        size_t i = 0;
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
            //deplrec.compute(m_msh, cl);
            //stab.compute(m_msh, cl, deplrec.oper);

            hyperelasticity.compute(m_msh, cl, lf, g, boundary_neumann, gradrec.oper, m_solution_data.at(i), m_elas_param);
            //hyperelasticity.compute2(m_msh, cl, lf, g, boundary_neumann, deplrec.oper, m_solution_data.at(i), m_elas_param);

            if(m_elas_param.adaptative_stab){
               m_beta = m_adpat_stab[i];
            }
            else
               m_beta = m_elas_param.tau;
               
            /////// NON LINEAIRE /////////
            tc.tic();
            dynamic_matrix<scalar_type> lhs = hyperelasticity.K_int + m_beta * stab.data;

            dynamic_vector<scalar_type> rhs = hyperelasticity.RTF  -  m_beta * stab.data * m_solution_data.at(i) ;
            dynamic_vector<scalar_type> rhs_cell = rhs.block(0,0, cbs, 1);


            dynamic_vector<scalar_type> x = statcond.recover(m_msh, cl, lhs, rhs_cell, xFs);
            m_postprocess_data.push_back(x);

            i++;
        }
        tc.toc();

        pi.time_postprocess = tc.to_double();

        return pi;
    }



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
      if(iter == 0)
         initial_residual = m_system_rhs.norm();

      scalar_type residual = m_system_rhs.norm();
      scalar_type max_error = m_system_rhs.maxCoeff();

      scalar_type relative_error = residual / initial_residual;

      if(initial_residual == scalar_type{0.0})
      {
         relative_error = 0.0;
         max_error = 0.0;
      }

      if(m_verbose){
         std::string s_iter = "   " + std::to_string(iter) + "               ";
         s_iter.resize(9);

         if(iter == 0){
            std::cout << "--------------------------------------------------------------" << std::endl;
            std::cout << "| Iteration |  Residual l2  | Relative error | Maximum error |" << std::endl;
            std::cout << "--------------------------------------------------------------" << std::endl;

         }
         std::ios::fmtflags f( std::cout.flags() );
         std::cout.precision(5);
         std::cout.setf(std::iostream::scientific, std::iostream::floatfield);
         std::cout << "| " << s_iter << " |   " << residual << " |   " << relative_error << "  |  " << max_error << "  |" << std::endl;
         std::cout << "--------------------------------------------------------------" << std::endl;
         std::cout.flags( f );
      }

      if(relative_error <= epsilon || residual <= epsilon ){
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
      solution_cells.clear();
      solution_cells = m_solution_cells;
      assert(m_solution_cells.size() == solution_cells.size());

      solution_faces.clear();
      solution_faces = m_solution_faces;
      assert(m_solution_faces.size() == solution_faces.size());

      solution_lagr.clear();
      solution_lagr = m_solution_lagr;
      assert(m_solution_lagr.size() == solution_lagr.size());

      solution.clear();
      solution = m_solution_data;
      assert(m_solution_data.size() == solution.size());
    }
};
