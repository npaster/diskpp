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
#include "hho_nl_diffusion.hpp"
#include "NewtonSolver_diffusion/newton_solver_diffusion.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>


struct offline_info
{
    double  time_offline;
};

struct solve_info
{
    double  time_solver;
};


template<typename Mesh, typename Point>
class NL_diffusion_solver
{
   typedef Mesh                                       mesh_type;
   typedef typename mesh_type::scalar_type            scalar_type;
   typedef typename mesh_type::cell                   cell_type;
   typedef typename mesh_type::face                   face_type;

   typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
   typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;

   typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
   typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

   typedef dynamic_matrix<scalar_type>         matrix_dynamic;
   typedef dynamic_vector<scalar_type>         vector_dynamic;

   size_t m_cell_degree, m_face_degree, m_degree;

   const mesh_type& m_msh;

   std::vector<vector_dynamic>                    m_solution_data;
   std::vector<matrix_dynamic>                    m_data_offline;
   std::vector<vector_dynamic>         m_solution_cells, m_solution_faces, m_solution_lagr;

   bool m_verbose;

public:
   NL_diffusion_solver(const mesh_type& msh, size_t degree, int l = 0)
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

   offline_info
   compute_offline()
   {
      disk::gradient_reconstruction< mesh_type,
                                          cell_basis_type,
                                          cell_quadrature_type,
                                          face_basis_type,
                                          face_quadrature_type> gradrec(m_degree);
                                          
                                          
      disk::diffusion_like_stabilization< mesh_type,
                                          cell_basis_type,
                                          cell_quadrature_type,
                                          face_basis_type,
                                          face_quadrature_type> stab(m_degree);

      m_data_offline.reserve(m_msh.cells_size());

      offline_info ai;
      bzero(&ai, sizeof(ai));

      timecounter tc;

      tc.tic();
      for (auto& cl : m_msh)
      {
         gradrec.compute(m_msh, cl);
         stab.compute(m_msh, cl, gradrec.oper);
         dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
         m_data_offline.push_back(loc);
      }
      tc.toc();
      ai.time_offline += tc.to_double();
      
      assert(m_data_offline.size() == m_msh.cells_size());

      return ai;
   }

   //template<typename DeplFunction, typename StressFunction>
   void
   compute_initial_state()//const DeplFunction& df, const StressFunction& bcf)
   {
      m_solution_data.clear();
      m_solution_cells.clear();
      m_solution_faces.clear();
      m_solution_lagr.clear();
      
      m_solution_data.reserve(m_msh.cells_size());
      m_solution_cells.reserve(m_msh.cells_size());
      m_solution_faces.reserve(m_msh.faces_size());
      m_solution_lagr.reserve(m_msh.boundary_faces_size());

      cell_basis_type cell_basis = cell_basis_type(m_degree);
      face_basis_type face_basis = face_basis_type(m_degree);

      const size_t num_cell_dofs = cell_basis.size();
      const size_t num_face_dofs = face_basis.size();

      for(auto& cl : m_msh)
      {
         auto fcs = faces(m_msh, cl);
         const size_t num_faces = fcs.size();
         m_solution_data.push_back(vector_dynamic::Zero(num_cell_dofs + num_faces * num_face_dofs));
         m_solution_cells.push_back(vector_dynamic::Zero(num_cell_dofs));
      }
      
      for(size_t i = 0; i < m_msh.faces_size(); i++)
      {
         m_solution_faces.push_back(vector_dynamic::Zero(num_face_dofs));
      }
      
      for(size_t i = 0; i < m_msh.boundary_faces_size(); i++){
         m_solution_lagr.push_back(vector_dynamic::Zero(num_face_dofs));   
      }    
   }

   template<typename LoadFunction, typename BoundaryConditionFunction>
   solve_info
   compute(const LoadFunction& lf, const BoundaryConditionFunction& bcf,
      size_t n_time_step = 1)
   {

      solve_info ai;
      bzero(&ai, sizeof(ai));

      timecounter tc;

      NewtonRaphson_solver_diffusion<Mesh> newton_solver(m_msh, m_degree, 0);
      
      newton_solver.initialize(m_solution_cells, m_solution_faces,
                                 m_solution_lagr, m_solution_data);

      const scalar_type delta_t = 1.0/n_time_step;

      scalar_type time = 0.0;

      for (size_t n = 0; n < n_time_step; n++)
      {
         time += delta_t;
         tc.tic();
         if(m_verbose){
            std::cout << "----------------------------------------------" << std::endl;
            std::cout << "************** Time step " << n+1 << "/" << n_time_step << " *****************" << std::endl;
         }

         auto rlf = [&lf, &time](const Point& p) -> auto {
            return time*lf(p);
         };

         auto rbcf = [&bcf, &time](const Point& p) -> auto {
            return time*bcf(p);
         };

         auto newton_info = newton_solver.compute(rlf, bcf, m_data_offline);

         tc.toc();
         ai.time_solver += tc.to_double();

         if(m_verbose){
            std::cout << "** Time in this step " << tc.to_double() << " sec" << std::endl;
            std::cout << "**** Assembly time: " << newton_info.time_assembly << " sec" << std::endl;
            std::cout << "**** Solver time: " << newton_info.time_solve << " sec" << std::endl;
            std::cout << "**** Postprocess time: " << newton_info.time_post << " sec" << std::endl;
         }
      }

      
      newton_solver.save_solutions(m_solution_cells, m_solution_faces,
                                 m_solution_lagr, m_solution_data);
      
      return ai;
   }



    template<typename AnalyticalSolution>
    scalar_type
    compute_l2_error(const AnalyticalSolution& as)
    {
        scalar_type err_dof = 0.0;

        disk::projector<mesh_type, cell_basis_type, cell_quadrature_type,
                        face_basis_type, face_quadrature_type> projk(m_degree);

        size_t i = 0;

        for (auto& cl : m_msh)
        {
            auto x = m_solution_cells.at(i++);
            dynamic_vector<scalar_type> true_dof = projk.compute_cell(m_msh, cl, as);
            dynamic_vector<scalar_type> comp_dof = x.block(0,0,true_dof.size(), 1);
            dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
            err_dof += diff_dof.dot(projk.cell_mm * diff_dof);
        }

        return sqrt(err_dof);
    }
//
//     void
//     plot_solution(const std::string& filename)
//     {
//         std::ofstream ofs(filename);
//
//         size_t cell_i = 0;
//         for (auto& cl : m_msh)
//         {
//             auto x = m_postprocess_data.at(cell_i++);
//             auto qps = m_bqd.cell_quadrature.integrate(m_msh, cl);
//             for (auto& qp : qps)
//             {
//                 auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, qp.point());
//
//                 scalar_type pot = 0.0;
//                 for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_cell_degree).size(); i++)
//                     pot += phi[i] * x(i);
//
//                 auto tp = qp.point();
//                 for (size_t i = 0; i < mesh_type::dimension; i++)
//                     ofs << tp[i] << " ";
//                 ofs << pot << std::endl;
//             }
//         }
//
//         ofs.close();
//     }
};
