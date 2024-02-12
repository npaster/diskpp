/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

/* this file deals with the heat equation with a backward Euler method */
#include <iostream>
#include <regex>
#include <sstream>
#include <iomanip>

#include <unistd.h>

#include "loaders/loader.hpp"
#include "methods/hho"
#include "solvers/solver.hpp"
#include "output/silo.hpp"
#include "timecounter.h"
#include "colormanip.h"
#include "bases/bases.hpp"

using namespace disk;
using namespace Eigen;

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
void
export_to_silo(const Mesh<T, 2, Storage>& msh,
               const Matrix<T, Dynamic, 1>& data, int cycle = -1)
{
    disk::silo_database silo;

    if (cycle == -1)
        silo.create("heat.silo");
    else
    {
        std::stringstream ss;
        ss << "out_" << cycle << ".silo";
        silo.create(ss.str());
    }

    silo.add_mesh(msh, "mesh");

    silo.add_variable("mesh", "sol", data, disk::zonal_variable_t );
    silo.close();
}

template<typename Mesh>
bool
unsteady_laplacian_solver(const Mesh& msh, size_t degree, typename Mesh::coordinate_type dt,
                          typename Mesh::coordinate_type penalization)
{
    typedef typename Mesh::coordinate_type  scalar_type;
    typedef typename Mesh::point_type       point_type;

    hho_degree_info hdi(degree, degree, degree+1);


    auto num_cells = msh.cells_size();
    auto num_faces = msh.faces_size();


    SparseMatrix<scalar_type>           LHS;
    Matrix<scalar_type, Dynamic, 1>     RHS;
    std::vector<Triplet<scalar_type>>   triplets;

    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto system_size = cbs*num_cells + fbs*num_faces;

    LHS = SparseMatrix<scalar_type>(system_size, system_size);
    RHS = Matrix<scalar_type, Dynamic, 1>::Zero(system_size);
    Matrix<scalar_type, Dynamic, 1> Mu_prev = Matrix<scalar_type, Dynamic, 1>::Zero(system_size);

    timecounter tc;

    int K = 1; // index for the exact solution

    auto rhs_fun = [](const point_type& pt) -> auto { return 0.0; };
    auto bcs_fun = [](const point_type& pt) -> auto { return 0.0; };
    auto solution= [K](const double t, const point_type& pt) -> auto { return std::exp(-2*M_PI*M_PI*K*K*t) * std::sin(M_PI*K*pt.x()) * std::sin(M_PI*K*pt.y()); };
    auto init_fun = [K](const point_type& pt) -> scalar_type {
			return std::sin(M_PI*K*pt.x()) * std::sin(M_PI*K*pt.y());
		    };

    tc.tic();
    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto fcs    = faces(msh, cl);
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);
        Matrix<scalar_type, Dynamic, Dynamic> lhs = gr.second + stab;

        Matrix<scalar_type, Dynamic, 1> rhs = Matrix<scalar_type, Dynamic, 1>::Zero(lhs.cols());
        rhs.head(cbs) = make_rhs(msh, cl, cb, rhs_fun, 1);

        apply_dirichlet_via_nitsche(msh, cl, gr.first, lhs, rhs, hdi, bcs_fun, penalization);

        Matrix<scalar_type, Dynamic, Dynamic> cell_mass = make_mass_matrix(msh, cl, cb);
        lhs = lhs*dt;
        lhs.block(0,0,cbs,cbs) += cell_mass;

        /* Make local-to-global mapping */
        std::vector<size_t> l2g(cbs + fcs.size() * fbs);
        for (size_t i = 0; i < cbs; i++)
            l2g[i] = cell_i*cbs + i;

        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto f_ofs = offset(msh, fcs[i]);
            for (size_t j = 0; j < fbs; j++)
                l2g[cbs + i*fbs+j] = cbs*num_cells + f_ofs*fbs+j;
        }

        /* Assemble cell contributions */
        for (size_t i = 0; i < lhs.rows(); i++)
        {
            for (size_t j = 0; j < lhs.cols(); j++)
                triplets.push_back( Triplet<scalar_type>(l2g[i], l2g[j], lhs(i,j)) );

            RHS(l2g[i]) += rhs(i);
        }

	/* set initial datum */
        Mu_prev.block(cell_i*cbs,0,cbs,1) = cell_mass * disk::project_function(msh, cl, degree, init_fun);

        cell_i++;
    }

    LHS.setFromTriplets(triplets.begin(), triplets.end());

    tc.toc();
    std::cout << " Assembly time: " << tc << std::endl;


    Matrix<scalar_type, Dynamic, 1> u;
    scalar_type t = 0.0;
    scalar_type L2H1_error = 0.;

    size_t freq_exp = 100;

    // time loop
    for (size_t i = 0; t < 2.0; i++, t += dt)
    {
        if(i % freq_exp == 0)
            std::cout << "Step " << i << std::endl;
        disk::solvers::pardiso_params<scalar_type> pparams;
        Matrix<scalar_type, Dynamic, 1> Mupp = Mu_prev;
        Mupp.tail(msh.faces_size() * fbs) = Matrix<scalar_type, Dynamic, 1>::Zero(msh.faces_size() * fbs);
        Mupp = Mupp + dt*RHS;
        mkl_pardiso(pparams, LHS, Mupp, u);

        Matrix<scalar_type, Dynamic, 1> sol_silo = Matrix<scalar_type, Dynamic, 1>::Zero(msh.cells_size());

        for (size_t i = 0; i < msh.cells_size(); i++)
            sol_silo(i) = u(i*cbs);

        if(i % freq_exp == 0)
            export_to_silo( msh, sol_silo, i );

        Mu_prev = Matrix<scalar_type, Dynamic, 1>::Zero(system_size);

        cell_i = 0;
        for (auto& cl : msh)
        {
            auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            Matrix<scalar_type, Dynamic, Dynamic> cell_mass = make_mass_matrix(msh, cl, cb);
            Mu_prev.block(cbs*cell_i, 0, cbs, 1) = cell_mass*u.block(cbs*cell_i, 0, cbs, 1);


            // solution
            auto sol_fun = [solution,t](const point_type& pt) -> scalar_type { return solution(t,pt); };

            /* compute L2-H1-error of the current time step */
            Matrix<scalar_type, Dynamic, 1> sol_ex
                = disk::project_function(msh, cl, degree, sol_fun);
            Matrix<scalar_type, Dynamic, 1> diff = u.block(cell_i*cbs, 0, cbs, 1) - sol_ex;

            const auto qps = integrate(msh, cl, 2*hdi.cell_degree());
            for(auto& qp : qps)
            {
                const auto g_phi = cb.eval_gradients( qp.point() );
                Matrix<scalar_type, 1, 2> grad = Matrix<scalar_type, 1, 2>::Zero();
                for (size_t i = 0; i < cbs; i++ )
                    grad += diff(i) * g_phi.block(i, 0, 1, 2);

                L2H1_error += dt * qp.weight() * grad.dot(grad);
            }

            cell_i++;
        }
   }

    std::cout << "L2-H1-error = " << std::sqrt(L2H1_error) << std::endl;

    /*
    size_t systsz = LHS.rows();
    size_t nnz = LHS.nonZeros();


    Matrix<scalar_type, Dynamic, 1> sol = Matrix<scalar_type, Dynamic, 1>::Zero(systsz);

    tc.tic();
    disk::solvers::pardiso_params<scalar_type> pparams;
    mkl_pardiso(pparams, LHS, RHS, sol);
    tc.toc();
    std::cout << " Solution time: " << tc << std::endl;

    Matrix<scalar_type, Dynamic, 1> sol_silo = Matrix<scalar_type, Dynamic, 1>::Zero(msh.cells_size());

    for (size_t i = 0; i < msh.cells_size(); i++)
        sol_silo(i) = sol(i*cbs);

    export_to_silo( msh, sol_silo );
    */

    return true;
}

/* run main with :
   ./unsteady_laplacian -m ../../../diskpp/meshes/2D_quads/diskpp/testmesh-16-16.quad -k 1 -t 0.01
*/


int main(int argc, char **argv)
{
    using T = double;
    disk::cartesian_mesh<T, 2> msh;

    T           dt = 0.01;
    size_t      degree = 1;
    char *      mesh_filename = nullptr;
    int ch;
    while ( (ch = getopt(argc, argv, "k:m:t:")) != -1 )
    {
	switch(ch)
        {
            case 'k':
                degree = std::stoi(optarg);
                break;

            case 'm':
                mesh_filename = optarg;
                break;

            case 't':
                dt = std::stof(optarg);
                break;

            default:
                std::cout << "Invalid option" << std::endl;
                return 1;
	}
    }

    msh = load_cartesian_2d_mesh<T>(mesh_filename);

    std::cout << "Mesh loaded ..." << std::endl;
    unsteady_laplacian_solver(msh, degree, dt, 1.0);
}

