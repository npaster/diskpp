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
#include "hho/hho_vector_bq.hpp"
#include "hho/hho_nl_vector.hpp"
#include "hho/hho_bq.hpp"
#include "../ElasticityParameters.hpp"
#include "../BoundaryConditions.hpp"
#include "../Informations.hpp"
#include "../hyperelasticity2_elementary_computation.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>



template<typename BQData>
class NewtonRaphson_step_hyperelasticity2
{
   typedef typename BQData::mesh_type          mesh_type;
   typedef typename mesh_type::scalar_type     scalar_type;

   typedef dynamic_matrix<scalar_type>         matrix_dynamic;
   typedef dynamic_vector<scalar_type>         vector_dynamic;


   typedef disk::gradient_reconstruction_elas_full_bq<BQData>               gradrec_type;

   typedef Hyperelasticity::Hyperelasticity2<BQData>                         hyperelasticity_type;

   typedef disk::diffusion_like_static_condensation_bq<BQData>              statcond_type;
   typedef disk::assembler_nl_vector_bq<BQData>                             assembler_type;


   typedef typename assembler_type::sparse_matrix_type       sparse_matrix_type;
   typedef typename assembler_type::vector_type              vector_type;

    #ifdef HAVE_INTEL_MKL
        typedef Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver_type;
    #else
        typedef Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver_type;
    #endif

    sparse_matrix_type     m_system_matrix;
    vector_type            m_system_rhs, m_system_solution;

    solver_type             solver;

    const BQData&                               m_bqd;

    const mesh_type& m_msh;

    std::vector<vector_dynamic>        m_postprocess_data, m_solution_data;
    std::vector<vector_dynamic>         m_solution_cells, m_solution_faces, m_solution_lagr;

    bool m_verbose;

    BoundaryConditions m_boundary_condition;

    ElasticityParameters m_elas_param;

    scalar_type initial_residual;


public:
   NewtonRaphson_step_hyperelasticity2(const mesh_type& msh, const BQData& bqd, const ElasticityParameters elas_param,
                                       const BoundaryConditions& boundary_conditions)
   : m_msh(msh), m_verbose(false), m_elas_param(elas_param), m_bqd(bqd),
    m_boundary_condition(boundary_conditions)
    {}

    bool    verbose(void) const     { return m_verbose; }
    void    verbose(bool v)         { m_verbose = v; }



    void
    initialize( const std::vector<vector_dynamic>& initial_solution_cells,
                const std::vector<vector_dynamic>& initial_solution_faces,
                const std::vector<vector_dynamic>& initial_solution_lagr,
                const std::vector<vector_dynamic>& initial_solution)
    {
      m_solution_cells.clear();
      m_solution_cells = initial_solution_cells;
      assert(m_msh.cells_size() == m_solution_cells.size());

      m_solution_faces.clear();
      m_solution_faces = initial_solution_faces;
      assert(m_msh.faces_size() == m_solution_faces.size());

      m_solution_lagr.clear();
      m_solution_lagr = initial_solution_lagr;

      m_solution_data.clear();
      m_solution_data = initial_solution;
      assert(m_msh.cells_size() == m_solution_data.size());
    }

    template<typename LoadFunction, typename BoundaryConditionFunction, typename NeumannFunction>
    AssemblyInfo
    assemble(const LoadFunction& lf, const BoundaryConditionFunction& bcf, const NeumannFunction& g)
    {
        gradrec_type gradrec(m_bqd);
        hyperelasticity_type hyperelasticity(m_bqd);

        statcond_type statcond(m_bqd);

        assembler_type assembler(m_msh, m_bqd, m_boundary_condition);

        AssemblyInfo ai;

        timecounter tc, ttot;

        ttot.tic();
        size_t i = 0;

        for (auto& cl : m_msh)
        {
            tc.tic();
            gradrec.compute(m_msh, cl, false);
            tc.toc();
            ai.m_time_gradrec += tc.to_double();

            tc.tic();
            hyperelasticity.compute(m_msh, cl, lf, g, m_boundary_condition.boundary_neumann(), gradrec.oper(), m_solution_data.at(i), m_elas_param);


            /////// NON LINEAIRE /////////
            dynamic_matrix<scalar_type> lhs = hyperelasticity.K_int;
            dynamic_vector<scalar_type> rhs = hyperelasticity.RTF;

            tc.toc();
            ai.m_time_elem += tc.to_double();
            ai.m_time_law += hyperelasticity.time_law;
            tc.tic();
            auto scnp = statcond.compute(m_msh, cl, lhs, rhs, true);
            tc.toc();
            ai.m_time_statcond += tc.to_double();

            assembler.assemble(m_msh, cl, scnp);

            i++;
        }

         assembler.impose_boundary_conditions(m_msh, bcf, m_solution_faces, m_solution_lagr, m_boundary_condition);
         assembler.finalize(m_system_matrix, m_system_rhs);

         ttot.toc();
         ai.m_time_assembly = ttot.to_double();
         ai.m_linear_system_size = m_system_matrix.rows();
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

    size_t
    test_aurrichio(void)
    {

      Eigen::SelfAdjointEigenSolver<sparse_matrix_type> es;
      es.compute(m_system_matrix);

      if(es.info() != Eigen::Success) {
          std::cerr << "ERROR: Could not compute eigenvalues of the matrix" << std::endl;
      }

      auto ev = es.eigenvalues();

      size_t nb_negative_eigenvalue(0);

      size_t i_ev = 0;
      while(ev(i_ev) <= scalar_type(0.0))
      {
         nb_negative_eigenvalue++;
         i_ev++;
      }

      const scalar_type min_ev = ev.minCoeff();
      const scalar_type max_ev = ev.maxCoeff();

      if(m_verbose){
         std::cout << "******* Eigenvalues test ********" << std::endl;
         std::cout << "Number of eigenvalues: " << ev.size()  << std::endl;
         std::cout << "Number of negative eigenvalues: " << nb_negative_eigenvalue  << std::endl;
         std::cout << "Maximum eigenvalue: " << max_ev << std::endl;
         std::cout << "Minimum eigenvalue: " << min_ev << std::endl;
         std::cout << "Conditionning number: " << std::abs(max_ev/min_ev) << std::endl;
      }

      return nb_negative_eigenvalue;
    }

    SolveInfo
    solve(void)
    {
       #ifdef HAVE_INTEL_MKL
       Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
       #else
       Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
       #endif

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

       return SolveInfo(m_system_matrix.rows(), m_system_matrix.nonZeros(), tc.to_double());
    }



    template<typename LoadFunction, typename NeumannFunction>
    PostprocessInfo
    postprocess(const LoadFunction& lf, const NeumannFunction& g)
    {
        gradrec_type gradrec(m_bqd);
        hyperelasticity_type hyperelasticity(m_bqd);
        statcond_type statcond(m_bqd);

        const size_t fbs = (m_bqd.face_basis.range(0, m_bqd.face_degree())).size();
        const size_t cbs = (m_bqd.cell_basis.range(0, m_bqd.cell_degree())).size();

        PostprocessInfo pi;

        m_postprocess_data.clear();

        m_postprocess_data.reserve(m_msh.cells_size());

        timecounter tc, ttot;
        tc.tic(); ttot.tic();

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

            tc.tic();
            gradrec.compute(m_msh, cl, false);
            tc.toc();
            pi.m_time_gradrec += tc.to_double();

            tc.tic();

            hyperelasticity.compute(m_msh, cl, lf, g, m_boundary_condition.boundary_neumann(), gradrec.oper(), m_solution_data.at(i), m_elas_param);

            /////// NON LINEAIRE /////////
            dynamic_matrix<scalar_type> lhs = hyperelasticity.K_int;

            dynamic_vector<scalar_type> rhs = hyperelasticity.RTF;
            dynamic_vector<scalar_type> rhs_cell = rhs.block(0,0, cbs, 1);

            tc.toc();
            pi.m_time_elem += tc.to_double();
            pi.m_time_law += hyperelasticity.time_law;

            tc.tic();
            dynamic_vector<scalar_type> x = statcond.recover(m_msh, cl, lhs, rhs_cell, xFs);
            tc.toc();
            pi.m_time_statcond += tc.to_double();

            m_postprocess_data.push_back(x);

            i++;
        }
        ttot.toc();

        pi.m_time_post = ttot.to_double();

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

        size_t cbs = (m_bqd.cell_basis.range(0, m_bqd.cell_degree())).size();

        assert(m_postprocess_data.size() == m_solution_cells.size());

        for(size_t i=0; i < m_postprocess_data.size(); i++){
            assert(m_solution_cells.at(i).size() == cbs);
            m_solution_cells.at(i) += (m_postprocess_data.at(i)).block(0,0,cbs,1);
        }

        size_t fbs = m_bqd.face_basis.size();


        for(size_t i=0; i < m_solution_faces.size(); i++){
            assert(m_solution_faces.at(i).size() == fbs);
            m_solution_faces.at(i) += m_system_solution.block(i * fbs, 0, fbs, 1);
        }

        const size_t lagrange_offset = m_solution_faces.size() * fbs;
        const size_t num_lagr_dofs = fbs/m_msh.dimension;
        for(size_t i=0; i < m_solution_lagr.size(); i++){
           const size_t pos = lagrange_offset + m_boundary_condition.begin_lag_conditions_faceI(i);
           const size_t size = num_lagr_dofs * m_boundary_condition.nb_lag_conditions_faceI(i);
           m_solution_lagr.at(i) += m_system_solution.block(lagrange_offset + i * fbs, 0, fbs, 1);
        }
    }



    bool test_convergence(const scalar_type epsilon, const size_t iter, scalar_type& error)
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

       size_t nb_faces_dof = m_bqd.face_basis.size() * m_msh.faces_size();

       scalar_type norm_depl = (m_system_rhs.head(nb_faces_dof)).norm();
       scalar_type norm_lag = (m_system_rhs.tail(m_system_rhs.size() - nb_faces_dof)).norm();

       if(m_verbose){
          std::string s_iter = "   " + std::to_string(iter) + "               ";
          s_iter.resize(9);

          if(iter == 0){
             std::cout << "------------------------------------------------------------------------------------------------" << std::endl;
             std::cout << "| Iteration |  Residual l2  | Relative error | Maximum error | Residual face  |  Residual BC   |" << std::endl;
             std::cout << "------------------------------------------------------------------------------------------------" << std::endl;

          }
          std::ios::fmtflags f( std::cout.flags() );
          std::cout.precision(5);
          std::cout.setf(std::iostream::scientific, std::iostream::floatfield);
          std::cout << "| " << s_iter << " |   " << residual << " |   " << relative_error << "  |  " << max_error << "  |   "
          << norm_depl << "  |   " << norm_lag << "  |" << std::endl;
          std::cout << "------------------------------------------------------------------------------------------------" << std::endl;
          std::cout.flags( f );
       }

       error = std::min(max_error, residual);

       if(error <= epsilon ){
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
