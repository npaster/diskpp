#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <unistd.h>

#include <map>

#include "colormanip.h"

#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"
#include "methods/hho"
#include "solvers/solver.hpp"

#include "../tests/common.hpp"

/***************************************************************************/
/* RHS definition */
template<typename Mesh>
struct rhs_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 1, Storage> >
{
    typedef Mesh<T,1,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        return M_PI * M_PI * sin_px;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;

    scalar_type
    operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        return 2.0 * M_PI * M_PI * sin_px * sin_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct rhs_functor<Mesh<T, 3, Storage>>
{
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;

    scalar_type
    operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto sin_pz = std::sin(M_PI * pt.z());
        return 3.0 * M_PI * M_PI * sin_px * sin_py * sin_pz;
    }
};

template<typename Mesh>
auto
make_rhs_function(const Mesh& msh)
{
    return rhs_functor<Mesh>();
}

/***************************************************************************/
/* Expected solution definition */
template<typename Mesh>
struct solution_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor< Mesh<T, 1, Storage> >
{
    typedef Mesh<T,1,Storage>               mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    scalar_type operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        return sin_px;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor< Mesh<T, 2, Storage> >
{
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;

    scalar_type
    operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        return sin_px * sin_py;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct solution_functor<Mesh<T, 3, Storage>>
{
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;

    scalar_type
    operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto sin_pz = std::sin(M_PI * pt.z());
        return sin_px * sin_py * sin_pz;
    }
};

template<typename Mesh>
auto
make_solution_function(const Mesh& msh)
{
    return solution_functor<Mesh>();
}

/***************************************************************************/
/* Expected solution definition */
template<typename Mesh>
struct gradient_functor;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct gradient_functor<Mesh<T, 2, Storage>>
{
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;
    typedef disk::static_vector<scalar_type, 2> result_type;

    result_type
    operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto cos_px = std::cos(M_PI * pt.x());
        auto cos_py = std::cos(M_PI * pt.y());

        auto gx = cos_px * sin_py;
        auto gy = sin_px * cos_py;

        return M_PI * result_type{gx, gy};
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct gradient_functor<Mesh<T, 3, Storage>>
{
    typedef Mesh<T, 3, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;
    typedef disk::static_vector<scalar_type, 3> result_type;

    result_type
    operator()(const point_type& pt) const
    {
        auto sin_px = std::sin(M_PI * pt.x());
        auto sin_py = std::sin(M_PI * pt.y());
        auto sin_pz = std::sin(M_PI * pt.z());
        auto cos_px = std::cos(M_PI * pt.x());
        auto cos_py = std::cos(M_PI * pt.y());
        auto cos_pz = std::cos(M_PI * pt.z());

        auto gx = cos_px * sin_py * sin_pz;
        auto gy = sin_px * cos_py * sin_pz;
        auto gz = sin_px * sin_py * cos_pz;
        return M_PI * result_type{gx, gy, gz};
    }
};

template<typename Mesh>
auto
make_gradient_function(const Mesh& msh)
{
    return gradient_functor<Mesh>();
}

using namespace disk;

template<typename Mesh>
typename Mesh::coordinate_type
run_hho_diffusion_solver(const Mesh& msh, size_t degree, const StabSize stab_diam_F, bool print = false)
{
    using T = typename Mesh::coordinate_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              vector_type;
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_type;

    hho_degree_info hdi(degree, degree, degree);

    if (print)
    {
        std::cout << "Mesh has " << msh.cells_size() << " elements." << std::endl;
    }

    const size_t odi = 2;

    auto rhs_fun  = make_rhs_function(msh);
    auto sol_fun  = make_solution_function(msh);
    auto grad_fun = make_gradient_function(msh);

    auto assembler = make_diffusion_assembler(msh, hdi);

    for (auto& cl : msh)
    {
        auto        cb   = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto        gr   = make_scalar_hho_laplacian(msh, cl, hdi);
        auto        stab = make_scalar_hho_stabilization(msh, cl, gr.first, hdi, stab_diam_F);
        auto        rhs  = make_rhs(msh, cl, cb, rhs_fun, odi);
        matrix_type A    = gr.second + stab;
        auto        sc   = make_scalar_static_condensation(msh, cl, hdi, A, rhs);
        assembler.assemble(msh, cl, sc.first, sc.second, sol_fun, odi);
    }

    assembler.finalize();

    size_t systsz = assembler.LHS.rows();
    size_t nnz    = assembler.LHS.nonZeros();

    if (print)
    {
        std::cout << "System has " << assembler.LHS.rows() << " unknowns and ";
        std::cout << assembler.LHS.nonZeros() << " nonzeros." << std::endl;
    }

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(systsz);

    disk::solvers::pardiso_params<T> pparams;

    if (print)
    {
        std::cout << "Running pardiso" << std::endl;
        pparams.report_factorization_Mflops = true;
    }

    mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

    T           error      = 0.0;
    vector_type flux_faces = vector_type::Zero(msh.faces_size());

    for (auto& cl : msh)
    {
        auto        cb   = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto        gr   = make_scalar_hho_laplacian(msh, cl, hdi);
        auto        stab = make_scalar_hho_stabilization(msh, cl, gr.first, hdi, stab_diam_F);
        auto        rhs  = make_rhs(msh, cl, cb, rhs_fun, odi);
        matrix_type A    = gr.second + stab;

        vector_type locsol = assembler.take_local_data(msh, cl, sol, sol_fun, odi);

        vector_type fullsol = make_scalar_static_decondensation(msh, cl, hdi, A, rhs, locsol);

        vector_type realsol = project_function(msh, cl, hdi, sol_fun, odi);

        auto diff = realsol - fullsol;
        error += diff.dot(A * diff);

        // L2-norm
        // matrix_type MM = make_mass_matrix(msh, cl, cb);
        // error += diff.head(cb.size()).dot(MM * diff.head(cb.size()));

        // compute fluxes
        const vector_type grad_u = gr.first * fullsol;
        const auto        gb     = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
        const auto        gbs    = gb.size() - 1;

        assert(grad_u.size() == gbs);

        const auto fcs = faces(msh, cl);

        const auto        adjoint = make_scalar_hho_stabilization_adjoint(msh, cl, gr.first, hdi, stab_diam_F);
        const vector_type flux_u  = adjoint * fullsol;

        size_t fc_off = 0;

        for (auto& fc : fcs)
        {
            const auto no = normal(msh, cl, fc);

            const auto diff_deg = std::max(hdi.face_degree(), hdi.cell_degree());
            const auto db       = make_scalar_monomial_basis(msh, fc, diff_deg);
            const auto dbs      = db.size();

            const vector_type flux_uF = flux_u.segment(fc_off, dbs);

            // compute fluxes F = grad(u).n + coeff_stab*(uF-uT)
            T flux_grad = T(0);
            T flux_stab = T(0);

            const auto qpf = integrate(msh, fc, std::max(hdi.grad_degree(), diff_deg));

            for (auto& qp : qpf)
            {
                // eval velocity
                const Eigen::Matrix<T, Eigen::Dynamic, Mesh::dimension> gphi =
                  gb.eval_gradients(qp.point()).block(1, 0, gbs, Mesh::dimension);
                const auto grad = eval(grad_u, gphi);
                flux_grad += qp.weight() * grad.dot(no);

                // eval stabilization term
                const auto dphi = db.eval_functions(qp.point());
                flux_stab += qp.weight() * (eval(flux_uF, dphi));
            }

            const auto face_id = msh.lookup(fc);
            flux_faces(face_id) += flux_grad + flux_stab;
            fc_off += dbs;

            if (msh.is_boundary(face_id))
            {
                // T flux_diri = T(0);
                // // accuracy of fluxes is really dependent of the integration order
                // const auto qpg = integrate(msh, fc, hdi.grad_degree());

                // const auto grad_u_real = project_function(msh, fc, hdi.grad_degree(), grad_fun, odi);
                // const auto pb          = make_vector_monomial_basis(msh, fc, hdi.grad_degree());

                // for (auto& qp : qpg)
                // {
                //     const auto gphi = pb.eval_functions(qp.point());
                //     const auto grad = eval(grad_u_real, gphi);
                //     // flux_diri += qp.weight() * grad_fun(qp.point()).dot(-no);
                //     flux_diri += qp.weight() * grad.dot(-no);
                // }

                // flux_faces(face_id) += flux_diri;

                // For the moment I do not know how to compute Dirichlet boundary
                flux_faces(face_id) = 0.;
            }
        }
    }

    if (print)
    {
        std::cout << "h = " << disk::average_diameter(msh) << " ";
        std::cout << "err = " << std::sqrt(error) << std::endl;
        std::cout << "Equilibrated fluxes -> norm: " << flux_faces.norm() << ", sum: " << flux_faces.sum() << std::endl;
    }

    return std::sqrt(error);
}

template<typename Mesh>
struct test_functor
{
    /* Expect k+1 convergence (hho stabilization, energy norm) */
    typename Mesh::coordinate_type
    operator()(const Mesh& msh, size_t degree) const
    {
        return run_hho_diffusion_solver(msh, degree, StabSize::hR, false);
    }

    size_t
    expected_rate(size_t k)
    {
        return k + 1;
    }
};

int
main(int argc, char** argv)
{
    using T = double;

    size_t      num_elems = 16;
    size_t      degree = 1;
    char *      mesh_filename = nullptr;
    bool        stat_cond     = true, use_mesh;
    StabSize    stab_diam_F   = StabSize::hF;

    int ch;
    while ( (ch = getopt(argc, argv, "k:m:twN:")) != -1 )
    {
        switch (ch)
        {
            case 'c': stab_diam_F = StabSize::hT; break;

            case 'k': degree = std::stoi(optarg); break;

            case 'N':
                /* Only for 1D */
                num_elems = std::stoi(optarg);
                break;

            case 'm':
                mesh_filename = optarg;
                use_mesh      = true;
                break;

            case '?':
            default: std::cout << "Invalid option" << std::endl; return 1;
        }
    }

    if (use_mesh)
    {
        /* FVCA5 2D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$")))
        {
            std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
            disk::generic_mesh<T, 2> msh;
            disk::load_mesh_fvca5_2d(mesh_filename, msh);
            run_hho_diffusion_solver(msh, degree, stab_diam_F, true);
            return 0;
        }

        /* Netgen 2D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$")))
        {
            std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
            disk::simplicial_mesh<T, 2> msh;
            disk::load_mesh_netgen(mesh_filename, msh);
            run_hho_diffusion_solver(msh, degree, stab_diam_F, true);
            return 0;
        }

        /* DiSk++ cartesian 2D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.quad$")))
        {
            std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
            disk::cartesian_mesh<T, 2> msh;
            disk::load_mesh_diskpp_cartesian(mesh_filename, msh);
            run_hho_diffusion_solver(msh, degree, stab_diam_F, true);
            return 0;
        }

        /* Poly 2d*/
        if (std::regex_match(mesh_filename, std::regex(".*\\.poly2d$")))
        {
            std::cout << "Guessed mesh format: Poly2D format" << std::endl;
            disk::generic_mesh<double, 2> msh;
            disk::load_mesh_poly2d<double>(mesh_filename, msh);
            run_hho_diffusion_solver(msh, degree, stab_diam_F, true);
            return 0;
        }

        /* Netgen 3D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$")))
        {
            std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
            disk::simplicial_mesh<T, 3> msh;
            disk::load_mesh_netgen(mesh_filename, msh);
            run_hho_diffusion_solver(msh, degree, stab_diam_F, true);
            return 0;
        }

        /* DiSk++ cartesian 3D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.hex$")))
        {
            std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
            disk::cartesian_mesh<T, 3> msh;
            disk::load_mesh_diskpp_cartesian(mesh_filename, msh);
            run_hho_diffusion_solver(msh, degree, stab_diam_F, true);
            return 0;
        }

        /* FVCA6 3D */
        if (std::regex_match(mesh_filename, std::regex(".*\\.msh$")))
        {
            std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
            disk::generic_mesh<T, 3> msh;
            disk::load_mesh_fvca6_3d<T>(mesh_filename, msh);
            run_hho_diffusion_solver(msh, degree, stab_diam_F, true);
            return 0;
        }
    }
    else
    {
        tester<test_functor> tstr;
        tstr.run();
        return 0;
    }
}