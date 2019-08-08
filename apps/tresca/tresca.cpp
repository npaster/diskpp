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

template<typename Mesh>
void
renumber_boundaries_gv2d(Mesh& msh)
{
    using T      = typename Mesh::coordinate_type;
    auto storage = msh.backend_storage();

    size_t b1 = 0;
    size_t b2 = 0;
    size_t b3 = 0;
    size_t b4 = 0;
    size_t b5 = 0;
    size_t b6 = 0;
    size_t b7 = 0;

    for (size_t face_i = 0; face_i < msh.faces_size(); face_i++)
    {
        auto fc = *std::next(msh.faces_begin(), face_i);
        if (storage->boundary_info.at(face_i).is_boundary)
        {
            const auto bar = barycenter(msh, fc);
            if (bar.x() >= 16.68 && bar.x() <= 37.3)
            {
                if (bar.y() >= T(8.07) && bar.y() <= T(8.10))
                {
                    storage->boundary_info.at(face_i).boundary_id = 1;
                    b4++;
                }
                else if (bar.y() <= T(5.30) && bar.y() >= T(4.755))
                {
                    storage->boundary_info.at(face_i).boundary_id = 2;
                    b5++;
                }
                else
                {
                    storage->boundary_info.at(face_i).boundary_id = 7;
                    b7++;
                }
            }
            else if (bar.x() >= 37.3 && bar.x() <= 59.8)
            {
                if (bar.y() >= T(8.07) && bar.y() <= T(8.10))
                {
                    storage->boundary_info.at(face_i).boundary_id = 1;
                    b4++;
                }
                else if (bar.y() <= T(6.80) && bar.y() >= T(5.30))
                {
                    storage->boundary_info.at(face_i).boundary_id = 2;
                    b5++;
                }
                else
                {
                    storage->boundary_info.at(face_i).boundary_id = 7;
                    b7++;
                }
            }
            else if (bar.x() >= 63.96 && abs(bar.y()) < T(1E-8))
            {
                storage->boundary_info.at(face_i).boundary_id = 3;
                b3++;
            }
            else if (std::abs(bar.x() - 10.8811) <= T(1E-5))
            {
                storage->boundary_info.at(face_i).boundary_id = 4;
                b4++;
            }
            else if (bar.y() > -5 && abs(bar.x() - 66.77) < T(1E-5))
            {
                storage->boundary_info.at(face_i).boundary_id = 5;
                b5++;
            }
            else
            {
                storage->boundary_info.at(face_i).boundary_id = 6;
                b6++;
            }
        }
    }

    std::cout << "b1: " << b1 << std::endl;
    std::cout << "b2: " << b2 << std::endl;
    std::cout << "b3: " << b3 << std::endl;
    std::cout << "b4: " << b4 << std::endl;
    std::cout << "b5: " << b5 << std::endl;
    std::cout << "b6: " << b6 << std::endl;
    std::cout << "b7: " << b7 << std::endl;
}

template<typename Mesh>
void
renumber_boundaries_gv3d(Mesh& msh)
{
    using T      = typename Mesh::coordinate_type;
    auto storage = msh.backend_storage();

    size_t b1 = 0;
    size_t b2 = 0;
    size_t b3 = 0;
    size_t b4 = 0;
    size_t b5 = 0;
    size_t b6 = 0;
    size_t b7 = 0;
    size_t b8 = 0;
    size_t b9 = 0;

    for (size_t face_i = 0; face_i < msh.faces_size(); face_i++)
    {
        auto fc = *std::next(msh.faces_begin(), face_i);
        if (storage->boundary_info.at(face_i).is_boundary)
        {
            const auto bar = barycenter(msh, fc);
            const auto r   = std::sqrt(bar.y() * bar.y() + bar.z() * bar.z());
            if (std::abs(bar.z()) <= T(1E-8))
            {
                storage->boundary_info.at(face_i).boundary_id = 1;
                b1++;
            }
            else if (std::abs(bar.y()) <= T(1E-8))
            {
                storage->boundary_info.at(face_i).boundary_id = 2;
                b2++;
            }
            else if (std::abs(bar.x() - 10.8811) <= T(1E-5))
            {
                storage->boundary_info.at(face_i).boundary_id = 3;
                b3++;
            }
            else if (std::abs(bar.x() - 66.77) <= T(1E-5))
            {
                storage->boundary_info.at(face_i).boundary_id = 7;
                b7++;
            }
            else if (r > T(8.70))
            {
                storage->boundary_info.at(face_i).boundary_id = 9;
                b9++;
            }
            else if (bar.x() >= 16.70 && bar.x() <= 37.3)
            {
                if (r >= T(8.07) && r <= T(8.10))
                {
                    storage->boundary_info.at(face_i).boundary_id = 4;
                    b4++;
                }
                else if (r <= T(5.30) && r >= T(4.755))
                {
                    storage->boundary_info.at(face_i).boundary_id = 5;
                    b5++;
                }
                else
                {
                    b8++;
                }
            }
            else if (bar.x() >= 37.3 && bar.x() <= 59.8)
            {
                if (r >= T(8.07) && r <= T(8.10))
                {
                    storage->boundary_info.at(face_i).boundary_id = 4;
                    b4++;
                }
                else if (r <= T(6.80) && r >= T(5.30))
                {
                    storage->boundary_info.at(face_i).boundary_id = 5;
                    b5++;
                }
                else
                {
                    b8++;
                }
            }
            else
            {
                // std::cout << bar << std::endl;
                b6++;
                // throw std::invalid_argument("dont find the boundaries");
            }
        }
    }
    std::cout << "b1: " << b1 << std::endl;
    std::cout << "b2: " << b2 << std::endl;
    std::cout << "b3: " << b3 << std::endl;
    std::cout << "b4: " << b4 << std::endl;
    std::cout << "b5: " << b5 << std::endl;
    std::cout << "b6: " << b6 << std::endl;
    std::cout << "b7: " << b7 << std::endl;
    std::cout << "b8: " << b8 << std::endl;
    std::cout << "b9: " << b9 << std::endl;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void run_tresca_solver(Mesh<T, 2, Storage>& msh, const ParamRun<T>& rp, const disk::MaterialData<T>& material_data)
{
    typedef Mesh<T, 2, Storage>                         mesh_type;
    typedef static_vector<T, 2>                         result_type;
    typedef static_matrix<T, 2, 2>                      result_grad_type;
    typedef disk::vector_boundary_conditions<mesh_type> Bnd_type;

    auto load = [material_data](const point<T, 2>& p) -> result_type { return result_type{0, 0}; };

    auto solution = [material_data](const point<T, 2>& p) -> result_type { return result_type{0, 0}; };

    Bnd_type bnd(msh);
    // bnd.addDirichletEverywhere(solution);

    auto s = [rp](const point<T, 2>& p) -> T { return rp.m_threshold; };

    // Bostan2
    // auto zero = [material_data](const point<T, 2>& p) -> result_type { return result_type{0.0, 0}; };
    // auto neum = [material_data](const point<T, 2>& p) -> result_type { return result_type{400.0, 0}; };

    // bnd.addContactBC(disk::SIGNORINI_FACE, 6, s);
    // bnd.addDirichletBC(disk::DIRICHLET, 10, zero);
    // // bnd.addDirichletBC(disk::DX, 6, zero);
    // bnd.addNeumannBC(disk::NEUMANN, 3, neum);

    // // Hertz
    // auto depl = [material_data](const point<T, 2>& p) -> result_type { return result_type{0.0, -5.0}; };
    // auto zero = [material_data](const point<T, 2>& p) -> result_type { return result_type{0.0, 0.0}; };

    // bnd.addContactBC(disk::SIGNORINI_FACE, 26, s);
    // bnd.addContactBC(disk::SIGNORINI_FACE, 16, s);

    // bnd.addDirichletBC(disk::DX, 7, zero);
    // bnd.addDirichletBC(disk::DX, 14, zero);
    // bnd.addDirichletBC(disk::DY, 4, depl);
    // bnd.addDirichletBC(disk::DY, 19, depl);

    // GV
    renumber_boundaries_gv2d(msh);

    auto zero = [material_data](const point<T, 2>& p) -> result_type { return result_type{0.0, 0}; };
    auto neum = [material_data](const point<T, 2>& p) -> result_type {
        const result_type n     = result_type{0.0, 1.0};
        const result_type t     = result_type{1.0, 0};
        const T           coeff = 11200;

        if (p.x() >= 18.0 && p.x() <= 30.0)
        {
            return 1.2 * coeff * (-n);
        }
        else if (p.x() <= 50.0)
        {
            return coeff * (-n - 0.06 * t);
        }
        else if (p.x() <= 60.0)
        {
            return 0.45 * coeff * (-n - 0.06 * t);
        }

        return result_type::Zero();
    };

    bnd.addContactBC(disk::SIGNORINI_FACE, 1, s);
    bnd.addNeumannBC(disk::NEUMANN, 2, neum);
    bnd.addDirichletBC(disk::CLAMPED, 5, zero);
    bnd.addDirichletBC(disk::DY, 3, zero);

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
    nl.compute_equivalent_plastic_strain_GP("plastic2D_p_GP_cont.msh");
    nl.compute_is_plastic_GP("isplastic2D_GP_cont.msh");
    nl.compute_continuous_deformed("deformed2D_cont.msh");
    // nl.compute_discontinuous_deformed("deformed2D_disc.msh");
    nl.compute_continuous_pres("pres2D_cont.msh");
    nl.contact_quantities("contact.dat");
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
void run_tresca_solver(Mesh<T, 3, Storage>& msh, const ParamRun<T>& rp, const disk::MaterialData<T>& material_data)
{
    typedef Mesh<T, 3, Storage>                         mesh_type;
    typedef static_vector<T, 3>                         result_type;
    typedef static_matrix<T, 3, 3>                      result_grad_type;
    typedef disk::vector_boundary_conditions<mesh_type> Bnd_type;

    auto load = [material_data](const point<T, 3>& p) -> result_type { return result_type{0, 0, 0}; };

    auto solution = [material_data](const point<T, 3>& p) -> result_type { return result_type{0, 0, 0}; };

    Bnd_type bnd(msh);
    // bnd.addDirichletEverywhere(solution);

    auto s = [rp](const point<T, 3>& p) -> T { return rp.m_threshold; };

    // // Bostan2
    // auto zero = [material_data](const point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, 0.0}; };
    // auto neum = [material_data](const point<T, 3>& p) -> result_type { return result_type{400.0, 0.0, 0.0}; };

    // bnd.addContactBC(disk::SIGNORINI_FACE, 23);
    // bnd.addDirichletBC(disk::CLAMPED, 27, zero);
    // bnd.addDirichletBC(disk::DZ, 31, zero);
    // bnd.addDirichletBC(disk::DZ, 33, zero);
    // bnd.addNeumannBC(disk::NEUMANN, 3, neum);

    // // Hertz
    // auto depl= [material_data](const point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, -2.0}; };
    // auto zero = [material_data](const point<T, 3>& p) -> result_type { return result_type{0.0, 0.0, 0.0}; };

    // bnd.addContactBC(disk::SIGNORINI_FACE, 31, s);
    // bnd.addDirichletBC(disk::DX, 19, zero);
    // bnd.addDirichletBC(disk::DX, 37, zero);
    // bnd.addDirichletBC(disk::DY, 27, zero);
    // bnd.addDirichletBC(disk::DY, 40, zero);
    // bnd.addDirichletBC(disk::DZ, 14, depl);

    // Gv
    renumber_boundaries_gv3d(msh);
    auto zero = [material_data](const point<T, 3>& p) -> result_type { return result_type::Zero(); };
    auto depl = [material_data](const point<T, 3>& p) -> result_type { return result_type{-0.2, -1, -1}; };

    auto neum = [material_data](const point<T, 3>& p) -> result_type {
        const result_type vec    = result_type{0.0, p.y(), p.z()};
        const result_type normal = -vec / vec.norm();
        const result_type vx     = result_type{1.0, 0, 0};
        const T           coeff  = 12400;

        if (p.x() >= 22.0 && p.x() <= 36.0)
        {
            return 1.1 * coeff * (-normal - 0.02 * vx);
        }
        else if (p.x() <= 50.0)
        {
            return coeff * (-normal - 0.08 * vx);
        }
        else if (p.x() <= 60.0)
        {
            return 0.4 * coeff * (-normal - 0.08 * vx);
        }

        return result_type::Zero();
    };

    bnd.addDirichletBC(disk::DZ, 1, zero);
    bnd.addDirichletBC(disk::DY, 2, zero);
    // bnd.addDirichletBC(disk::DIRICHLET, 5, neum);
    bnd.addNeumannBC(disk::NEUMANN, 5, neum);
    bnd.addDirichletBC(disk::CLAMPED, 7, zero);
    bnd.addDirichletBC(disk::CLAMPED, 9, zero);
    bnd.addContactBC(disk::SIGNORINI_FACE, 4, s);

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

    // nl.compute_discontinuous_displacement("depl3D_disc.msh");
    nl.compute_continuous_displacement("depl3D_cont.msh");
    // nl.compute_discontinuous_stress("stress3D_disc.msh");
    nl.compute_continuous_stress("stress3D_cont.msh");
    nl.compute_stress_GP("stress3D_GP.msh");
    //    nl.compute_equivalent_plastic_strain_GP("plastic_p_GP_cont.msh");
    //    nl.compute_is_plastic_GP("isplastic_GP_cont.msh");
    nl.compute_continuous_deformed("deformed3D_cont.msh");
    nl.compute_continuous_pres("pres3D_cont.msh");
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

    // // Hertz elasticity
    // RealType E  = 70;
    // RealType nu = 0.3;

    // material_data.setMu(E, nu);
    // material_data.setLambda(E, nu);

    // GV
    RealType E  = 210E3;
    RealType nu = 0.3;

    material_data.setMu(E, nu);
    material_data.setLambda(E, nu);

    // Hertz plasticity
    // RealType E  = 70;
    // RealType nu = 0.3;

    // material_data.setMu(E, nu);
    // material_data.setLambda(E, nu);

    // material_data.setK(0.0);
    // material_data.setH(7.0);
    // material_data.setSigma_y0(2.5);

    // material_data.setMu(1.0);
    // material_data.setLambda(1000.0);

    // Bostan
    // material_data.setMu(1000, 0.3);
    // material_data.setLambda(1000, 0.3);

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
