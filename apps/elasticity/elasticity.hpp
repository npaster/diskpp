/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
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

#pragma once

#include "common/eigen.hpp"
#include "quadratures/quadratures.hpp"
#include "bases/bases.hpp"

#include "hho/hho.hpp"

namespace disk {

    template<typename CellBasisType, typename CellQuadType,
    typename FaceBasisType, typename FaceQuadType, typename Mesh,
         typename Function>
dynamic_vector<typename Mesh::scalar_type>
compute_rhs_whole(const Mesh& msh, const typename Mesh::cell& cl,
            const Function& f, size_t degree)
{
    typedef dynamic_vector<typename Mesh::scalar_type> vector_type;

    auto cell_basis     = CellBasisType(degree);
    auto cell_quad      = CellQuadType(2*degree);
    auto face_basis     = FaceBasisType(degree);
    auto face_quad      = FaceQuadType(2*degree);

    auto fcs = faces(msh, cl);

    vector_type ret = vector_type::Zero(cell_basis.size() + fcs.size() * face_basis.size());

    auto cell_quadpoints = cell_quad.integrate(msh, cl);
    for (auto& qp : cell_quadpoints)
    {
        auto phi = cell_basis.eval_functions(msh, cl, qp.point());
        auto fval = f(qp.point());
        for (size_t i = 0; i < cell_basis.size(); i++)
            ret(i) += qp.weight() * mm_prod(fval, phi.at(i));
    }

    size_t offset_face = cell_basis.size();

     for (auto& fc : fcs)
        {
            auto face_quadpoints = face_quad.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
                auto phi = face_basis.eval_functions(msh, fc, qp.point());
                auto fval = f(qp.point());
                for (size_t i = 0; i < face_basis.size(); i++)
                    ret(offset_face + i) += qp.weight() * mm_prod(fval, phi.at(i));
            }

            offset_face += face_basis.size();
        }

    return ret;
}


template<typename CellBasisType, typename CellQuadType,
typename Mesh,
typename Function>
dynamic_vector<typename Mesh::scalar_type>
compute_rhs_cell(const Mesh& msh, const typename Mesh::cell& cl,
                  const Function& f, size_t degree)
{
   typedef dynamic_vector<typename Mesh::scalar_type> vector_type;

   auto cell_basis     = CellBasisType(degree);
   auto cell_quad      = CellQuadType(2*degree);


   vector_type ret = vector_type::Zero(cell_basis.size());

   auto cell_quadpoints = cell_quad.integrate(msh, cl);
   for (auto& qp : cell_quadpoints)
   {
      auto phi = cell_basis.eval_functions(msh, cl, qp.point());
      auto fval = f(qp.point());
      for (size_t i = 0; i < cell_basis.size(); i++)
         ret(i) += qp.weight() * mm_prod(fval, phi.at(i));
   }

   return ret;
}


   template<typename Mesh, typename CellBasisType, typename CellQuadType,
    typename FaceBasisType, typename FaceQuadType, typename Solution>
class l2_error
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef CellBasisType                       cell_basis_type;
    typedef CellQuadType                        cell_quadrature_type;
    typedef FaceBasisType                       face_basis_type;
    typedef FaceQuadType                        face_quadrature_type;


    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;




    cell_basis_type                             cell_basis;
    cell_quadrature_type                        cell_quadrature;

    face_basis_type                             face_basis;
    face_quadrature_type                        face_quadrature;



    size_t                                      m_degree;

public:
    matrix_type     oper;
    matrix_type     data;

    l2_error()
        : m_degree(1)
    {
        cell_basis          = cell_basis_type(m_degree);
        cell_quadrature     = cell_quadrature_type(m_degree);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    l2_error(size_t degree)
        : m_degree(degree)
    {
        cell_basis          = cell_basis_type(m_degree);
        cell_quadrature     = cell_quadrature_type(m_degree);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }



    scalar_type compute(const mesh_type& msh, const cell_type& cl,
                 const vector_type& x, Solution sol)
    {

         disk::projector_elas<mesh_type,
                                                   cell_basis_type,
                                                   cell_quadrature_type,
                                                   face_basis_type,
                                                   face_quadrature_type> proj(m_degree);


            dynamic_vector<scalar_type> true_dof = proj.compute_cell(msh, cl, sol);
            dynamic_vector<scalar_type> comp_dof = x.block(0,0,true_dof.size(), 1);

//             std::cout << "error" << std::endl;
//             for (size_t i = 0; i < comp_dof.size(); i++) {
//                std::cout << true_dof(i) <<" vs "<< comp_dof(i)  << '\n';
//             }
            dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
            scalar_type err_dof = diff_dof.dot(proj.cell_mm * diff_dof);


            return err_dof;

    }
};

template<typename MeshType, typename Solution>
void
test_condensation(const MeshType& msh, const Solution& solution, const size_t degree)
{
   typedef MeshType                            mesh_type;

   typedef typename mesh_type::scalar_type     scalar_type;

   typedef typename mesh_type::cell            cell_type;
   typedef typename mesh_type::face            face_type;

   typedef disk::quadrature<mesh_type, cell_type>   cell_quadrature_type;
   typedef disk::quadrature<mesh_type, face_type>   face_quadrature_type;
   typedef disk::scaled_monomial_vector_basis<mesh_type, cell_type>  cell_basis_type;
   typedef disk::scaled_monomial_vector_basis<mesh_type, face_type>  face_basis_type;


   disk::projector_elas<mesh_type,
                        cell_basis_type,
                        cell_quadrature_type,
                        face_basis_type,
                        face_quadrature_type> proj(degree);



   disk::diffusion_like_static_condensation<mesh_type,
                                                   cell_basis_type,
                                                   cell_quadrature_type,
                                                   face_basis_type,
                                                   face_quadrature_type> statcond(degree);

    disk::assembler_elas<mesh_type, face_basis_type,
                          face_quadrature_type> assembler(msh, degree);


    disk::l2_error<mesh_type, cell_basis_type, cell_quadrature_type, face_basis_type,
                                                   face_quadrature_type, Solution> error(degree);


   const size_t dim = msh.dimension;

    timecounter tc;

    tc.tic();
    for (auto& cl : msh)
    {

        dynamic_vector<scalar_type> tmp = proj.compute_whole(msh, cl, solution);
        dynamic_matrix<scalar_type> loc = proj.whole_mm;

        dynamic_matrix<scalar_type> rhs = compute_rhs_whole<cell_basis_type, cell_quadrature_type,
                                                       face_basis_type, face_quadrature_type>(msh, cl, solution, degree);

        auto sc = statcond.compute(msh, cl, loc, rhs, true);

        assembler.assemble(msh, cl, sc);
    }

    assembler.impose_boundary_conditions(msh, solution);
    assembler.finalize();
    tc.toc();

    std::cout << "Assembly time: " << tc << " seconds." << std::endl;

    tc.tic();

    #ifdef HAVE_INTEL_MKL
    Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
#else
    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
#endif

    size_t systsz = assembler.matrix.rows();
    size_t nnz = assembler.matrix.nonZeros();

     std::cout << "Starting linear solver..." << std::endl;
     std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
     std::cout << " * Matrix fill: " << 100.0*double(nnz)/(systsz*systsz) << "%" << std::endl;


    solver.analyzePattern(assembler.matrix);
    solver.factorize(assembler.matrix);
    dynamic_vector<scalar_type> X = solver.solve(assembler.rhs);
//#endif

    tc.toc();
    std::cout << "Solver time: " << tc << " seconds." << std::endl;

    face_basis_type face_basis(degree);
    size_t fbs = face_basis.size();
    cell_basis_type cell_basis(degree);
    size_t cbs = cell_basis.size();

    scalar_type l2 = 0.0;
    size_t cont_cell = 0;

    for (auto& cl : msh)
    {
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(num_faces*fbs);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            dynamic_vector<scalar_type> xF = dynamic_vector<scalar_type>::Zero(fbs);
            xF = X.block(face_id * fbs, 0, fbs, 1);
            xFs.block(face_i * fbs, 0, fbs, 1) = xF;
        }

        dynamic_matrix<scalar_type> rhs = compute_rhs_whole<cell_basis_type, cell_quadrature_type,
                                                       face_basis_type, face_quadrature_type>(msh, cl, solution, degree);


//         for(size_t i = 0; i <   num_faces*fbs; i++)
//             std::cout << xFs(i) << " vs " << rhs(cbs+i) << std::endl;

        dynamic_vector<scalar_type> tmp = proj.compute_whole(msh, cl, solution);
        dynamic_matrix<scalar_type> loc = proj.whole_mm;


        dynamic_matrix<scalar_type> cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, solution, degree);


        dynamic_vector<scalar_type> x = statcond.recover(msh, cl, loc, cell_rhs, xFs);

        l2 += error.compute(msh, cl, x, solution);
        cont_cell += 1;

    }

    std::cout << "test condensation l2_error: " << sqrt(l2) << std::endl;

}




template<typename MeshType, typename Function, typename Solution>
void test_gradient(MeshType& msh, const Function& load, const Solution& solution, size_t degree)
{
    typedef MeshType                            mesh_type;

    typedef typename mesh_type::scalar_type     scalar_type;

    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef disk::quadrature<mesh_type, cell_type>   cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>   face_quadrature_type;
    typedef disk::scaled_monomial_vector_basis<mesh_type, cell_type>  cell_basis_type;
    typedef disk::scaled_monomial_vector_basis<mesh_type, face_type>  face_basis_type;

    if (degree < 1)
    {
        std::cout << "Only K > 0 for this problem" << std::endl;
        return;
    }


    disk::sgradient_reconstruction_elas<mesh_type,
                                        cell_basis_type,
                                        cell_quadrature_type,
                                        face_basis_type,
                                        face_quadrature_type> gradrec(degree);



    disk::elas_like_stabilization<mesh_type,
                                             cell_basis_type,
                                             cell_quadrature_type,
                                             face_basis_type,
                                             face_quadrature_type> stab(degree);

    disk::diffusion_like_static_condensation<mesh_type,
                                                   cell_basis_type,
                                                   cell_quadrature_type,
                                                   face_basis_type,
                                                   face_quadrature_type> statcond(degree);

    disk::assembler_elas<mesh_type, face_basis_type,
                          face_quadrature_type> assembler(msh, degree);


    disk::l2_error<mesh_type, cell_basis_type, cell_quadrature_type, face_basis_type,
                                                   face_quadrature_type, Solution> error(degree);

      disk::projector_elas<mesh_type,
                                                                                             cell_basis_type,
                                                                                             cell_quadrature_type,
                                                                                             face_basis_type,
                                                                                             face_quadrature_type> proj(degree);

    scalar_type mu      = 1.0;
    scalar_type lambda  = 0.0;
    const size_t DIM = msh.dimension;

    timecounter tc;

    tc.tic();
    for (auto& cl : msh)
    {
        gradrec.compute(msh, cl);
        stab.compute(msh, cl, gradrec.oper);

        auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, load, degree);

//         std::cout << "gradrec" << std::endl;
//         std::cout << gradrec.oper << std::endl;
//         std::cout << "stab" << std::endl;
//         std::cout << stab.data << std::endl;
//         std::cout << "rhs" << std::endl;
//         std::cout << cell_rhs << std::endl;
        //std::cout << cell_rhs << '\n';
        assert(cell_rhs.size() == DIM * binomial(degree + DIM, degree));
        dynamic_matrix<scalar_type> loc =  2.0* mu * gradrec.data +
                                            2.0*mu * stab.data;
        auto sc = statcond.compute(msh, cl, loc, cell_rhs);

//         std::cout << "mat_cond" << std::endl;
//         std::cout << sc.first << std::endl;
//         std::cout << "rhs_cond" << std::endl;
//         std::cout << sc.second << std::endl;


        assembler.assemble(msh, cl, sc);
    }

    assembler.impose_boundary_conditions(msh, solution);
    assembler.finalize();
    tc.toc();

    //std::cout << "Assembly time: " << tc << " seconds." << std::endl;

    tc.tic();

//#ifdef HAVE_SOLVER_WRAPPERS
//    agmg_solver<scalar_type> solver;
//    dynamic_vector<scalar_type> X = solver.solve(assembler.matrix, assembler.rhs);
//#else

#ifdef HAVE_INTEL_MKL
    Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
#else
    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
#endif

    size_t systsz = assembler.matrix.rows();
    size_t nnz = assembler.matrix.nonZeros();

     /*std::cout << "Starting linear solver..." << std::endl;
     std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
     std::cout << " * Matrix fill: " << 100.0*double(nnz)/(systsz*systsz) << "%" << std::endl;*/


    solver.analyzePattern(assembler.matrix);
    solver.factorize(assembler.matrix);
    dynamic_vector<scalar_type> X = solver.solve(assembler.rhs);
//#endif

    tc.toc();
    //std::cout << "Solver time: " << tc << " seconds." << std::endl;

    face_basis_type face_basis(degree);
    auto fbs = face_basis.size();

    scalar_type l2 = 0.0;
    size_t cont_cell = 0;

    for (auto& cl : msh)
    {
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(num_faces*fbs);



        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            dynamic_vector<scalar_type> xF = dynamic_vector<scalar_type>::Zero(fbs);
            xF = X.block(face_id * fbs, 0, fbs, 1);
            xFs.block(face_i * fbs, 0, fbs, 1) = xF;

//             std::cout << "sol" << std::endl;
//             std::cout << xFs << std::endl;
        }

        gradrec.compute(msh, cl);
        stab.compute(msh, cl, gradrec.oper);

        auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, load, degree);
        dynamic_matrix<scalar_type> loc = 2.0* mu * gradrec.data +
                                          2.0* mu * stab.data;

        dynamic_vector<scalar_type> x = statcond.recover(msh, cl, loc, cell_rhs, xFs);

        auto dof_true = proj.compute_whole(msh, cl, solution);

        std::cout << "STU" << std::endl;
        std::cout << stab.data * x << std::endl;
        std::cout << "USTU " << x.transpose() * stab.data * x << std::endl;

        std::cout << "STU true" << std::endl;
        std::cout << stab.data * dof_true << std::endl;
        std::cout << "USTU " << dof_true.transpose() * stab.data * dof_true << std::endl;

//         std::cout << "recover" << std::endl;
//         std::cout << x << std::endl;

        l2 += error.compute(msh, cl, x, solution);
        cont_cell += 1;

    }

    std::cout << "grad test l2_error: " << l2 << std::endl;

}

} // namespace disk
