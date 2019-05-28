/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
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

#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <unistd.h>

#include <map>

#include "Informations.hpp"
#include "Parameters.hpp"
#include "boundary_conditions/boundary_conditions.hpp"
#include "loaders/loader.hpp"
#include "mechanics/behaviors/laws/materialData.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "tresca_solver.hpp"

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_tresca_solver(const Mesh<T, 2, Storage>& msh, const ParamRun<T>& rp, const disk::MaterialData<T>& material_data)
{
    typedef Mesh<T, 2, Storage>                         mesh_type;
    typedef static_vector<T, 2>                         result_type;
    typedef static_matrix<T, 2, 2>                      result_grad_type;
    typedef disk::vector_boundary_conditions<mesh_type> Bnd_type;

    auto load = [material_data](const point<T, 2>& p) -> result_type { return result_type{0, 0}; };

    auto solution = [material_data](const point<T, 2>& p) -> result_type { return result_type{0, 0}; };

    Bnd_type bnd(msh);
    // bnd.addDirichletEverywhere(solution);

    // Bostan2
    // auto zero = [material_data](const point<T, 2>& p) -> result_type { return result_type{0.0, 0}; };
    // auto depl = [material_data](const point<T, 2>& p) -> result_type { return result_type{0.0, -0.2}; };
    // auto neum = [material_data](const point<T, 2>& p) -> result_type { return result_type{400.0, 0}; };

    // bnd.addContactBC(disk::SIGNORINI_CELL, 6);
    // bnd.addDirichletBC(disk::DIRICHLET, 10, depl);

    // // Hertz
    // auto depl = [material_data](const point<T, 2>& p) -> result_type { return result_type{0.0, -2.0}; };

    // bnd.addContactBC(disk::SIGNORINI_FACE, 6);
    // bnd.addDirichletBC(disk::DIRICHLET, 3, depl);

    // Conv
    auto depl = [material_data](const point<T, 2>& p) -> result_type {
        T ux = -0.2 * p.y() * p.y() * p.y() * (p.x() - 0.5);
        T uy = -0.05 * (1.0 + p.y()) * (p.x() - 0.5) * (p.x() - 0.5) - 0.01 * p.y();
        return result_type{ux, uy};
    };

    auto load2 = [material_data](const point<T, 2>& p) -> result_type {
        T x  = p.x();
        T y  = p.y();
        T y2 = y * y;

        T fx = material_data.getLambda() * (0.1 * x - 0.05) +
               2.0 * material_data.getMu() * (0.05 * x + 2 * y * (0.3 * x - 0.15) - 0.0025);
        T fy = 0.6 * material_data.getLambda() * y2 + 2.0 * material_data.getMu() * (0.3 * y2 + 0.05 * y + 0.05);
        return result_type{fx, fy};
    };

    auto s = [rp](const point<T, 2>& p) -> T {

        return rp.m_threshold;
    };

    bnd.addContactBC(disk::SIGNORINI_FACE, 2, s);
    bnd.addDirichletBC(disk::DIRICHLET, 1, depl);

    // solver

    tresca_solver<mesh_type> nl(msh, bnd, rp, material_data);

    if (nl.verbose())
    {
        std::cout << "Solving the problem ..." << '\n';
    }

    SolverInfo solve_info = nl.compute(load);

    if (nl.verbose())
    {
        solve_info.printInfo();
    }

    // nl.compute_discontinuous_displacement("depl2D_disc.msh");
    nl.compute_continuous_displacement("depl2D_cont.msh");
    // nl.compute_discontinuous_stress("stress2D_disc.msh");
    nl.compute_continuous_stress("stress2D_cont.msh");
    nl.compute_stress_GP("stress2D_GP.msh");
    nl.compute_continuous_deformed("deformed2D_cont.msh");
    // nl.compute_discontinuous_deformed("deformed2D_disc.msh");
    nl.contact_quantities("contact.dat");
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void
run_tresca_solver(const Mesh<T, 3, Storage>& msh, const ParamRun<T>& rp, const disk::MaterialData<T>& material_data)
{
    typedef Mesh<T, 3, Storage>                         mesh_type;
    typedef static_vector<T, 3>                         result_type;
    typedef static_matrix<T, 3, 3>                      result_grad_type;
    typedef disk::vector_boundary_conditions<mesh_type> Bnd_type;

    auto load = [material_data](const point<T, 3>& p) -> result_type { return result_type{0, 0, 0}; };

    auto solution = [material_data](const point<T, 3>& p) -> result_type { return result_type{0, 0, 0}; };

    Bnd_type bnd(msh);
    // bnd.addDirichletEverywhere(solution);

    // // Bostan2
    // auto zero = [material_data](const point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, 0.0}; };
    // auto neum = [material_data](const point<T, 3>& p) -> result_type { return result_type{400.0, 0.0, 0.0}; };

    // bnd.addContactBC(disk::SIGNORINI_CELL, 23);
    // bnd.addDirichletBC(disk::CLAMPED, 27, zero);
    // bnd.addDirichletBC(disk::DZ, 31, zero);
    // bnd.addDirichletBC(disk::DZ, 33, zero);
    // bnd.addNeumannBC(disk::NEUMANN, 3, neum);

    auto s = [rp](const point<T, 3>& p) -> T {

        return rp.m_threshold;
    };

    // Hertz
    auto depl= [material_data](const point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, -5.0}; };

    bnd.addContactBC(disk::SIGNORINI_CELL, 17, s);
    bnd.addDirichletBC(disk::DIRICHLET, 13, depl);

    // Solver

    tresca_solver<mesh_type> nl(msh, bnd, rp, material_data);

    if (nl.verbose())
    {
        std::cout << "Solving the problem ..." << '\n';
    }

    SolverInfo solve_info = nl.compute(load);

    if (nl.verbose())
    {
        solve_info.printInfo();
    }

    nl.compute_discontinuous_displacement("depl3D_disc.msh");
    nl.compute_continuous_displacement("depl3D_cont.msh");
    // nl.compute_discontinuous_stress("stress3D_disc.msh");
    nl.compute_continuous_stress("stress3D_cont.msh");
    nl.compute_stress_GP("stress3D_GP.msh");
    // nl.compute_continuous_deformed("deformed3D_cont.msh");
    // nl.compute_discontinuous_deformed("deformed3D_disc.msh");
    nl.contact_quantities("contact.dat");
}

int
main(int argc, char** argv)
{
    using RealType = double;

    char* mesh_filename = nullptr;

    ParamRun<RealType> rp;

    const RealType MPa = 10E6;
    const RealType GPa = 10E9;

    // Elasticity Parameters
    disk::MaterialData<RealType> material_data;

    RealType E  = 1000;
    RealType nu = 0.4999;

    material_data.setMu(E, nu);
    material_data.setLambda(E, nu);
    // material_data.setMu(1.0);
    // material_data.setLambda(1000.0);

    int ch;

    while ((ch = getopt(argc, argv, "r:")) != -1)
    {
        switch (ch)
        {
            case 'r':
                if (!rp.readParameters(optarg))
                    exit(1);
                break;
            default: std::cout << "wrong arguments" << std::endl; exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    if (argc == 0)
    {
        std::cout << "Error" << std::endl;
        return 0;
    }

    mesh_filename = argv[0];

    /* Medit 2d*/
    if (std::regex_match(mesh_filename, std::regex(".*\\.medit2d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        auto msh = disk::load_medit_2d_mesh<RealType>(mesh_filename);
        run_tresca_solver(msh, rp, material_data);
        return 0;
    }

    /* Netgen 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$")))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
        auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
        run_tresca_solver(msh, rp, material_data);
        return 0;
    }

    /* Medit 3d*/
    if (std::regex_match(mesh_filename, std::regex(".*\\.medit3d$")))
    {
        std::cout << "Guessed mesh format: Medit format" << std::endl;
        auto msh = disk::load_medit_3d_mesh<RealType>(mesh_filename);
        run_tresca_solver(msh, rp, material_data);
        return 0;
    }
}
