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
#include "BehaviorLaws/behavior_laws.hpp"
#include "BehaviorLaws/material_parameters.hpp"
//#include "contrib/sol2/sol.hpp"

//#define USE_BLAS
#define FILL_COLMAJOR

namespace disk {

template<typename Mesh, typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class gradient_reconstruction_elas
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



    typedef material_tensor<scalar_type, mesh_type::dimension, mesh_type::dimension>
                                                material_tensor_type;

    cell_basis_type                             cell_basis;
    cell_quadrature_type                        cell_quadrature;

    face_basis_type                             face_basis;
    face_quadrature_type                        face_quadrature;

    size_t                                      m_degree;

public:
    matrix_type     oper;
    matrix_type     data;

    gradient_reconstruction_elas()
        : m_degree(1)
    {
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    gradient_reconstruction_elas(size_t degree)
        : m_degree(degree)
    {
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    void compute(const mesh_type& msh, const cell_type& cl)
    {
        material_tensor_type id_tens;
        id_tens = material_tensor_type::Identity();
        compute(msh, cl, id_tens);
    }

    void compute(const mesh_type& msh, const cell_type& cl,
                 const material_tensor_type& mtens)
    {
       const size_t DIM= msh.dimension;
       const size_t dpk1 = DIM * binomial(m_degree +1 + DIM, m_degree +1);
       const size_t dpk0 = DIM * binomial( DIM, 0);
       const size_t dpk = DIM * binomial(m_degree  + DIM, m_degree );
       const size_t dpkf = DIM * binomial(m_degree  + DIM -1, m_degree );

        matrix_type stiff_mat = matrix_type::Zero(cell_basis.size(), cell_basis.size());

        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
        for (auto& qp : cell_quadpoints)
        {
            auto dphi = cell_basis.eval_gradients(msh, cl, qp.point());
            assert(cell_basis.size() == dphi.size());
            assert(dpk1 == dphi.size());
            for(size_t i=0; i< dphi.size(); i++){
               for(size_t j=0; j< dphi.size(); j++){
               //    std::cout << dphi.at(i) << '\n';
               //    std::cout << dphi.at(j) << '\n';
               //   std::cout << i << " " << j << " " <<  mm_prod(dphi.at(i), dphi.at(j)) << '\n';
                  stiff_mat(i,j) += qp.weight() * mm_prod(dphi.at(i), dphi.at(j));
               }
            }
            // std::cout << " " << '\n';
        }

        /* LHS: take basis functions derivatives from degree 1 to K+1 */
        auto MG_rowcol_range = cell_basis.range(1, m_degree+1);
        assert(MG_rowcol_range.size() == (dpk1 - dpk0));
        matrix_type MG = take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

        /* RHS, volumetric part. */
        auto BG_row_range = cell_basis.range(1, m_degree+1);
        auto BG_col_range = cell_basis.range(0, m_degree);

        assert(BG_row_range.size() == (dpk1 - dpk0));
        assert(BG_col_range.size() == dpk) ;

        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        auto num_cell_dofs = cell_basis.range(0, m_degree).size();
        auto num_face_dofs = face_basis.size();

        assert(num_cell_dofs == dpk);
        assert(num_face_dofs == dpkf);

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        assert(dsr.total_size() == (num_cell_dofs + num_faces * num_face_dofs));

        matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

        BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
                                take(stiff_mat, BG_row_range, BG_col_range);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto current_face_range = dsr.face_range(face_i);
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);
            auto face_quadpoints = face_quadrature.integrate(msh, fc);

            for (auto& qp : face_quadpoints)
            {
                auto c_phi =
                    cell_basis.eval_functions(msh, cl, qp.point()); // 0, m_degree);
                auto c_dphi =
                    cell_basis.eval_gradients(msh, cl, qp.point()); // 1, m_degree+1);

                decltype(c_phi) c_dphi_n;

                assert(BG_row_range.from() == dpk0);
                assert(BG_row_range.to() == dpk1);

                for(size_t i=BG_row_range.from(); i< BG_row_range.to(); i++){
                  c_dphi_n.push_back(mm_prod(c_dphi.at(i) , n));
                }

                assert(c_dphi_n.size() == (dpk1 - dpk0));

                matrix_type  T= matrix_type::Zero(BG.rows(), BG_col_range.size());

                assert(c_dphi_n.size() == BG.rows());

                for(size_t i=0; i< BG.rows(); i++){
                  for(size_t j=0; j<BG_col_range.size(); j++){
                     T(i,j) = qp.weight() * mm_prod(c_dphi_n.at(i), c_phi.at(j));
                  }
               }


                BG.block(0, 0, BG.rows(), BG_col_range.size()) -= T;

                auto f_phi = face_basis.eval_functions(msh, fc, qp.point());

                assert(f_phi.size() == dpkf);

                matrix_type  F= matrix_type::Zero(BG.rows(), current_face_range.size());

                for(size_t i=0; i< BG.rows(); i++){
                  for(size_t j=0; j<current_face_range.size(); j++){
                     F(i,j) = qp.weight() * mm_prod(c_dphi_n.at(i), f_phi.at(j));
                  }
               }


                BG.block(0, current_face_range.min(),
                         BG.rows(), current_face_range.size()) += F;
            }
        }

        assert(MG.rows() ==MG.cols());
        assert(MG.cols() == BG.rows());

        oper  = MG.llt().solve(BG);    // GT
        data  = BG.transpose() * oper;  // A
    }
};



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
        auto fcs = faces(msh, cl);
        std::vector<size_t> l2g(fcs.size() * face_basis.size());
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


template<typename CellBasisType, typename Mesh>
static_matrix<typename Mesh::scalar_type,3,3>
compute_gradient_pt(const Mesh& msh, const typename Mesh::cell& cl,
            const dynamic_vector<typename Mesh::scalar_type>& gradrec_coeff,
            const point<typename Mesh::scalar_type,3>& pt, const size_t degree)
{
    typedef typename Mesh::scalar_type               scalar_type;
    typedef static_matrix<scalar_type,3,3>            gradient_value_type;

    CellBasisType cell_basis          = CellBasisType(degree+1);

    auto G_range = cell_basis.range(1, degree+1);

    gradient_value_type ret = gradient_value_type::Zero(3,3);

    auto dphi = cell_basis.eval_gradients(msh, cl, pt);
    assert(cell_basis.size() == dphi.size());
    assert((G_range.to()-G_range.from()) == gradrec_coeff.size());

    for(size_t i = G_range.from(); i< G_range.to(); i++){
        ret += gradrec_coeff(i - G_range.from()) * dphi.at(i);
    }

    return ret;
}

template<typename CellBasisType, typename Mesh>
static_matrix<typename Mesh::scalar_type,2,2>
compute_gradient_pt(const Mesh& msh, const typename Mesh::cell& cl,
            const dynamic_vector<typename Mesh::scalar_type>& gradrec_coeff,
            const point<typename Mesh::scalar_type,2>& pt, const size_t degree)
{
    typedef typename Mesh::scalar_type               scalar_type;
    typedef static_matrix<scalar_type,2,2>            gradient_value_type;

    CellBasisType cell_basis          = CellBasisType(degree+1);

    auto G_range = cell_basis.range(1, degree+1);

    gradient_value_type ret = gradient_value_type::Zero(2,2);

    auto dphi = cell_basis.eval_gradients(msh, cl, pt);
    assert(cell_basis.size() == dphi.size());
    assert((G_range.to()-G_range.from()) == gradrec_coeff.size());

    for(size_t i = G_range.from(); i< G_range.to(); i++){
        ret += gradrec_coeff(i - G_range.from())  * dphi.at(i);
    }

    return ret;
}

// The gradient of transformation tensor
template<typename T>
static_matrix<T,3,3>
compute_fgradient_pt(const static_matrix<T,3,3>& gradient)
{
    return gradient + static_matrix<T,3,3>::Identity();
}

template<typename T>
static_matrix<T,2,2>
compute_fgradient_pt(const static_matrix<T,2,2>& gradient)
{
    return gradient + static_matrix<T,2,2>::Identity();
}

// The right Cauch-Green tensor
template<typename T>
static_matrix<T,3,3>
compute_cgright_pt(const static_matrix<T,3,3>& fgradient)
{
    return fgradient.transpose() * fgradient;
}

template<typename T>
static_matrix<T,2,2>
compute_cgright_pt(const static_matrix<T,2,2>& fgradient)
{
    return fgradient.transpose() * fgradient;
}


template<typename CellBasisType, typename CellQuadType, typename Mesh>
std::pair<dynamic_matrix<typename Mesh::scalar_type>, dynamic_vector<typename Mesh::scalar_type> >
compute_elem(const Mesh& msh, const typename Mesh::cell& cl,
            const dynamic_vector<typename Mesh::scalar_type>& gradrec_coeff, const size_t degree)
{
    typedef typename Mesh::scalar_type  scalar_type;
    typedef dynamic_vector<scalar_type> vector_type;
    typedef dynamic_matrix<scalar_type> matrix_type;

    CellBasisType cell_basis          = CellBasisType(degree+1);
    CellQuadType cell_quadrature      = CellQuadType(2*(degree+1));

    auto G_range = cell_basis.range(1, degree+1);

    size_t dim_mat = G_range.to()-G_range.from();

    assert(gradrec_coeff.size() == dim_mat);

    matrix_type K = matrix_type::Zero(dim_mat,dim_mat);
    vector_type R = vector_type::Zero(dim_mat);

    //  a automatiser
    StVenantKirchhoffLaw<scalar_type>  law(1.0,1.0);

    auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
    for (auto& qp : cell_quadpoints)
    {
        auto dphi = cell_basis.eval_gradients(msh, cl, qp.point());
        //Compoute G(u)
        auto gradu = compute_gradient_pt<CellBasisType, Mesh>(msh, cl, gradrec_coeff, qp.point(), degree);
        //compute F(u) and F(u)^T
        auto fu = compute_fgradient_pt<scalar_type>(gradu);
        auto transpo_fu = fu.transpose();

        //compute transpo_fu * dphi
        decltype(dphi) fdphi;
        fdphi.reserve(dphi.size());

        for(size_t i = 0; i < dphi.size(); i++)
            fdphi.push_back(transpo_fu * dphi.at(i));

        //calculer C et S
        auto cu = compute_cgright_pt<scalar_type>(fu);

        auto pk2 = law.compute_PK2(cu);

        //compute C(u) : transpo_fu * dphi
        decltype(dphi) cdphi;
        cdphi.reserve(dphi.size());

        for(size_t i = 0; i < fdphi.size(); i++)
            cdphi.push_back( /* ajouter C */  fdphi.at(i));

        //compure K_geom
        for (size_t i = G_range.from(); i < G_range.to(); i++)
            for (size_t j = i; j < G_range.to(); j++)
                K(i - G_range.from(), j - G_range.from()) += qp.weight() * mm_prod(fdphi.at(i), cdphi.at(j));


        //compute  dphi * S(u)
        decltype(dphi) sdphi;
        sdphi.reserve(dphi.size());

        for(size_t i = 0; i < dphi.size(); i++)
            sdphi.push_back(dphi.at(i) * pk2);

        //compure K_elas
        for (size_t i = G_range.from(); i < G_range.to(); i++)
            for (size_t j = i; j < G_range.to(); j++)
                K(i - G_range.from(), j - G_range.from()) += qp.weight() * mm_prod(sdphi.at(i), dphi.at(j));

        //compute F(u) * S(u)
        auto fsu = fu /* * pk2 */ ;

        //compure R_int
        for (size_t i = G_range.from(); i < G_range.to(); i++)
            R(i - G_range.from()) += qp.weight() * mm_prod(fsu, dphi.at(i));
    }

    //we use the symetrie of K
    for (size_t i = 0; i < dim_mat; i++)
         for (size_t j = i + 1; j < dim_mat; j++)
             K(j,i) = K(i,j);


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

    return (GT.transpose() * Rint + coeff_stab * STu - Rext);
}


} // namespace disk
