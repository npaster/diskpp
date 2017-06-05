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
#include "bases/bases_utils.hpp"
#include "bases/bases_ranges.hpp"
#include "timecounter.h"
#include "hho/hho.hpp"
#include "BehaviorLaws/behaviorlaws.hpp"
#include "BehaviorLaws/maths_tensor.hpp"
#include "BehaviorLaws/material_parameters.hpp"
//#include "contrib/sol2/sol.hpp"

//#define USE_BLAS
#define FILL_COLMAJOR

namespace disk {

    template<typename Mesh, //typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class assembler_nl_elas
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    //typedef CellBasisType                       cell_basis_type;
    //typedef CellQuadType                        cell_quadrature_type;
    typedef FaceBasisType                       face_basis_type;
    typedef FaceQuadType                        face_quadrature_type;

    typedef dynamic_matrix<scalar_type>         matrix_type;


    //cell_basis_type                             cell_basis;
    //cell_quadrature_type                        cell_quadrature;

    face_basis_type                             face_basis;
    face_quadrature_type                        face_quadrature;

    size_t                                      m_degree;

    typedef Eigen::Triplet<scalar_type>         triplet_type;

    std::vector<triplet_type>                   m_triplets;
    size_t                                      m_num_unknowns;

public:

    typedef Eigen::SparseMatrix<scalar_type>    sparse_matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

    sparse_matrix_type      matrix;
    vector_type             rhs;

    assembler_nl_elas()                 = delete;

    assembler_nl_elas(const mesh_type& msh, size_t degree)
        : m_degree(degree)
    {
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);

        m_num_unknowns = face_basis.size() * (msh.faces_size() + msh.boundary_faces_size());
        matrix = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
        rhs = vector_type::Zero(m_num_unknowns);
    }

    template<typename LocalContrib>
    void
    assemble(const mesh_type& msh, const cell_type& cl, const LocalContrib& lc)
    {

       const size_t DIM= msh.dimension;
       const size_t dpk1 = DIM * binomial(m_degree +1 + DIM, m_degree +1);
       const size_t dpk0 = DIM * binomial( DIM, 0);
       const size_t dpk = DIM * binomial(m_degree  + DIM, m_degree );
       const size_t dpkf = DIM * binomial(m_degree  + DIM -1, m_degree );


       auto fcs = faces(msh, cl);
       std::vector<size_t> l2g(fcs.size() * face_basis.size());

       assert(dpkf == face_basis.size());

       for (size_t face_i = 0; face_i < fcs.size(); face_i++)
       {
          auto fc = fcs[face_i];
          auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
          if (!eid.first)
             throw std::invalid_argument("This is a bug: face not found");

          auto face_id = eid.second;

          auto face_offset = face_id * face_basis.size();

          auto pos = face_i * face_basis.size();

          for (size_t i = 0; i < face_basis.size(); i++)
             l2g[pos+i] = face_offset+i;
       }

       assert(lc.first.rows() == fcs.size()*dpkf);
       assert(lc.first.rows() == lc.first.cols());
       assert(lc.first.rows() == lc.second.size());
       assert(lc.second.size() == l2g.size());

       //std::cout << lc.second.size() << " " << l2g.size() << std::endl;

       #ifdef FILL_COLMAJOR
       for (size_t j = 0; j < lc.first.cols(); j++)
       {
          for (size_t i = 0; i < lc.first.rows(); i++)
             m_triplets.push_back( triplet_type( l2g.at(i), l2g.at(j), lc.first(i,j) ) );

          rhs(l2g.at(j)) += lc.second(j);
       }
       #else
       for (size_t i = 0; i < lc.first.rows(); i++)
       {
          for (size_t j = 0; j < lc.first.cols(); j++)
             m_triplets.push_back( triplet_type( l2g.at(i), l2g.at(j), lc.first(i,j) ) );

          rhs(l2g.at(i)) += lc.second(i);
       }
       #endif
    }

    template<typename Function>
    void
    impose_boundary_conditions(const mesh_type& msh, const Function& bc, const std::vector<vector_type>& sol_faces,
                               const std::vector<vector_type>& sol_lagr)
    {
        size_t fbs = face_basis.size();
        size_t face_i = 0;
        size_t lagr_i = 0;
        for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
        {
            auto bfc = *itor;

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            auto face_offset = face_id * fbs;
            auto face_offset_lagrange = (msh.faces_size() + face_i) * fbs;

            auto fqd = face_quadrature.integrate(msh, bfc);

            matrix_type MFF     = matrix_type::Zero(fbs, fbs);
            vector_type rhs_f   = vector_type::Zero(fbs);

            for (auto& qp : fqd)
            {
                auto f_phi = face_basis.eval_functions(msh, bfc, qp.point());
                for(size_t i=0; i< fbs; i++)
                   for(size_t j=0; j< fbs; j++)
                        MFF(i,j) += qp.weight() * mm_prod(f_phi.at(i), f_phi.at(j));

                for(size_t i=0; i< fbs; i++)
                  rhs_f(i) += qp.weight() * mm_prod( f_phi.at(i),  bc(qp.point()));

            }

            rhs_f -= MFF * sol_faces.at(face_id);

            vector_type rhs_l = MFF * sol_lagr.at(lagr_i);

#ifdef FILL_COLMAJOR
            for (size_t j = 0; j < MFF.cols(); j++)
            {
                for (size_t i = 0; i < MFF.rows(); i++)
                {
                    m_triplets.push_back( triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i,j)) );
                    m_triplets.push_back( triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i,j)) );
                }
                rhs(face_offset_lagrange+j) = rhs_f(j);
                rhs(face_offset+j) -= rhs_l(j);
            }
#else
            for (size_t i = 0; i < MFF.rows(); i++)
            {
                for (size_t j = 0; j < MFF.cols(); j++)
                {
                    m_triplets.push_back( triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i,j)) );
                    m_triplets.push_back( triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i,j)) );
                }
                rhs(face_offset_lagrange+i) = rhs_f(i);
                rhs(face_offset+j) -= rhs_l(j);
            }
#endif
            face_i++;
            lagr_i++;
        }
    }

    void
    finalize()
    {
        matrix.setFromTriplets(m_triplets.begin(), m_triplets.end());
        m_triplets.clear();
    }

    void
    finalize(sparse_matrix_type& mat, vector_type& vec)
    {
        mat = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
        mat.setFromTriplets(m_triplets.begin(), m_triplets.end());
        m_triplets.clear();
        vec = rhs;
    }
};


template<typename CellBasisType, typename CellQuadType, typename Mesh,
         typename Function>
dynamic_vector<typename Mesh::scalar_type>
compute_rhs_ext(const Mesh& msh, const typename Mesh::cell& cl,
            const Function& f, const size_t degree,
            const dynamic_vector<typename Mesh::scalar_type>& solution)
{
    typedef dynamic_vector<typename Mesh::scalar_type> vector_type;

    vector_type ret = vector_type::Zero(solution.size());

    vector_type rhs_ext = compute_rhs<CellBasisType, CellQuadType>(msh, cl, f, degree);

    ret.block(0,0, rhs_ext.rows(), rhs_ext.cols()) = rhs_ext;

    return ret;
}


template<typename GradBasisType, typename Mesh>
static_matrix<typename Mesh::scalar_type,3,3>
compute_gradient_pt(const Mesh& msh, const typename Mesh::cell& cl,
            const dynamic_vector<typename Mesh::scalar_type>& gradrec_coeff,
            const point<typename Mesh::scalar_type,3>& pt, const size_t degree)
{
    typedef typename Mesh::scalar_type               scalar_type;
    typedef static_matrix<scalar_type,3,3>            gradient_value_type;

    GradBasisType grad_basis          = GradBasisType(degree);

    const size_t DIM = msh.dimension;
    const size_t tpk = DIM * DIM * binomial(degree  + DIM, degree );

    gradient_value_type ret = gradient_value_type::Zero(3,3);

    auto gphi = grad_basis.eval_functions(msh, cl, pt);
    assert(gphi.size() == gradrec_coeff.size());
    assert(gphi.size() == tpk);

    for(size_t i = 0; i < gphi.size(); i++){
        ret += gradrec_coeff(i) * gphi.at(i);
    }

    return ret;
}


template<typename GradBasisType, typename Mesh>
static_matrix<typename Mesh::scalar_type,2,2>
compute_gradient_pt(const Mesh& msh, const typename Mesh::cell& cl,
                    const dynamic_vector<typename Mesh::scalar_type>& gradrec_coeff,
                    const point<typename Mesh::scalar_type,2>& pt, const size_t degree)
{
   typedef typename Mesh::scalar_type               scalar_type;
   typedef static_matrix<scalar_type,2,2>            gradient_value_type;
   
   GradBasisType grad_basis          = GradBasisType(degree);
   
   const size_t DIM = msh.dimension;
   const size_t tpk = DIM * DIM * binomial(degree  + DIM, degree );
   
   gradient_value_type ret = gradient_value_type::Zero(2,2);
   
   auto gphi = grad_basis.eval_functions(msh, cl, pt);
   assert(gphi.size() == gradrec_coeff.size());
   assert(gphi.size() == tpk);
   
   for(size_t i = 0; i < gphi.size(); i++){
      ret += gradrec_coeff(i) * gphi.at(i);
   }
   
   return ret;
}





template< typename GradBasisType, typename GradQuadType, typename Mesh>
std::pair<dynamic_matrix<typename Mesh::scalar_type>, dynamic_vector<typename Mesh::scalar_type> >
compute_elem(const Mesh& msh, const typename Mesh::cell& cl,
            const dynamic_vector<typename Mesh::scalar_type>& gradrec_coeff, const size_t degree)
{
    typedef typename Mesh::scalar_type  scalar_type;
    typedef dynamic_vector<scalar_type> vector_type;
    typedef dynamic_matrix<scalar_type> matrix_type;
    
    GradBasisType grad_basis          = GradBasisType(degree);
    GradQuadType grad_quadrature      = GradQuadType(2*degree);

    const size_t DIM= msh.dimension;
    const size_t dpk1 = DIM * binomial(degree +1 + DIM, degree +1);
    const size_t dpk0 = DIM * binomial( DIM, 0);
    const size_t dpk = DIM * binomial(degree  + DIM, degree );
    const size_t gpk = DIM * DIM * binomial(degree  + DIM, degree );
    const size_t dpkf = DIM * binomial(degree  + DIM -1, degree );

    auto G_range = grad_basis.range(0, degree);

    assert(G_range.from() == 0);
    assert(G_range.to() == dpk);
    assert(DIM *(G_range.size()) == gpk);

    size_t dim_mat = gpk;

    assert(gradrec_coeff.size() == dim_mat);

    matrix_type K = matrix_type::Zero(dim_mat,dim_mat);
    vector_type R = vector_type::Zero(dim_mat);

    //  a automatiser
    LaplaceLaw<scalar_type>  law(3);

    auto grad_quadpoints = grad_quadrature.integrate(msh, cl);
    for (auto& qp : grad_quadpoints)
    {
        auto gphi = grad_basis.eval_functions(msh, cl, qp.point());
        assert(gphi.size() == gpk);

        //Compoute G(u)
        auto gradu = compute_gradient_pt<GradBasisType, Mesh>(msh, cl, gradrec_coeff, qp.point(), degree);


//         std::cout << "Gu" << std::endl;
//         std::cout << gradu << std::endl;

        auto tensor_behavior = law.compute_whole(gradu);
        //std::cout << "PK2" << tensor_behavior.first << '\n';
        auto P = tensor_behavior.first;
        auto A = tensor_behavior.second;

//         std::cout << "P" << std::endl;
//         std::cout << P << std::endl;
//         std::cout << "A" << std::endl;
//         std::cout << A << std::endl;

//         auto eigen = Ce.eigenvalues();
//
//         std::cout << eigen << std::endl;

        //compute A(u) : gphi
        decltype(gphi) Agphi;
        Agphi.reserve(gphi.size());

        for(size_t i = 0; i < gphi.size(); i++)
        {
            Agphi.push_back(tm_prod(A, gphi.at(i)) );
        }

        //compure K
        for (size_t i = 0; i < gphi.size(); i++)
            for (size_t j = 0; j < gphi.size(); j++)
                K(i,j) += qp.weight() * mm_prod(gphi.at(i), Agphi.at(j));


        //compure R_int
        for (size_t i = 0; i < gphi.size(); i++)
            R(i) += qp.weight() * mm_prod(P, gphi.at(i));
    }

    return std::make_pair(K,R);
}



template<typename T>
dynamic_matrix<T>
assemble_lhs(const dynamic_matrix<T>& GT, const dynamic_matrix<T>& ST, const dynamic_matrix<T>& K,
              const T coeff_stab)
{
    assert(K.cols() == GT.rows());
    assert(K.rows() == K.cols());
    assert(GT.cols() == ST.rows());
    assert(ST.cols() == ST.rows());

    return (GT.transpose() * K * GT + coeff_stab * ST);
}

template<typename T>
dynamic_vector<T>
assemble_rhs(const dynamic_matrix<T>& GT, const dynamic_vector<T>& STu,
             const dynamic_vector<T>& Rint, const dynamic_vector<T>& Rext,
              const T coeff_stab)
{
    assert(GT.rows() == Rint.rows());
    assert(GT.cols() == STu.rows());
    assert(STu.rows() == Rext.rows());

//     std::cout << "Rint" << std::endl;
//     std::cout << Rint << std::endl;
//     std::cout << "Rext" << std::endl;
//     std::cout << Rext << std::endl;

    return (GT.transpose() * Rint + coeff_stab * STu - Rext);
}


template<typename Mesh, typename CellBasisType, typename CellQuadType,
typename FaceBasisType, typename FaceQuadType,
typename Solution>
void
test_gradient_reconstruction(const Mesh& msh, const typename Mesh::cell& cl,
                             const Solution& as, dynamic_matrix<typename Mesh::scalar_type> G,
                             size_t degree)
{
   disk::projector_elas<Mesh, CellBasisType, CellQuadType,
   FaceBasisType, FaceQuadType> projk(degree);

   disk::projector_elas<Mesh, CellBasisType, CellQuadType,
   FaceBasisType, FaceQuadType> projk1(degree+1);
   // projetion of the solution on P^(k)_d (solution of (u,v) = (f,v))
   auto proju = projk.compute_whole(msh, cl, as);

   //keep only the part of uT
   // multiply gradrec_oper(0 : N(k+1) -1, 0: Nvt) *  proju(1:N(k+1))
   auto Gu = (G * proju).eval();

   // projetion of the solution on P^(k+1)_d (solution of (u,v) = (f,v))
   auto sol = projk1.compute_cell(msh, cl, as);

   size_t dim = msh.dimension;



   double l2error = 0.0;

   for (size_t i = 0; i < Gu.size(); i++) {
      l2error +=  (Gu(i)-sol(i+dim))*( Gu(i)-sol(i+dim));
      //std::cout << sol(i+dim) <<  " vs " << Gu(i)  << " "  <<   Gu(i)-sol(i+dim)  <<  '\n';
   }

          std::cout <<  "l2error grad: "  <<   sqrt(l2error)  <<  '\n';

}


} // namespace disk
