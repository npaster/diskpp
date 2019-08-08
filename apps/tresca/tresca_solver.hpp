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

#include <iostream>
#include <sstream>

#include <list>
#include <vector>

#include "bases/bases.hpp"
#include "methods/hho"
#include "quadratures/quadratures.hpp"

#include "Informations.hpp"
#include "NewtonSolver/newton_solver.hpp"
#include "NewtonSolver/tresca_elementary_computation.hpp"
#include "Parameters.hpp"

#include "boundary_conditions/boundary_conditions.hpp"
#include "mechanics/behaviors/laws/behaviorlaws.hpp"
#include "mechanics/behaviors/laws/materialData.hpp"
#include "mechanics/behaviors/logarithmic_strain/LogarithmicStrain.hpp"
#include "mechanics/behaviors/tensor_conversion.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/gmsh_io.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

struct time_step
{
    double time;
    size_t level;
};

template<typename Mesh>
class tresca_solver
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef ParamRun<scalar_type>               param_type;
    typedef disk::MaterialData<scalar_type>     material_type;

    typedef dynamic_matrix<scalar_type> matrix_type;
    typedef dynamic_vector<scalar_type> vector_type;

    typedef disk::vector_boundary_conditions<mesh_type>           bnd_type;
    typedef disk::LinearIsotropicAndKinematicHardening<mesh_type> law_type;

    typename disk::hho_degree_info m_hdi;
    bnd_type                       m_bnd;
    const mesh_type&               m_msh;
    material_type                  m_material_data;

    law_type m_law;

    disk::gmsh_io<mesh_type> gmesh_io;

    std::vector<vector_type> m_solution, m_solution_faces;
    std::vector<matrix_type> m_gradient_precomputed, m_stab_precomputed;

    bool m_verbose, m_convergence;

    param_type m_rp;

    const static size_t dimension = mesh_type::dimension;

    int total_dof_depl_static;

    void
    init()
    {
        m_solution.clear();
        m_solution_faces.clear();

        m_solution.reserve(m_msh.cells_size());
        m_solution_faces.reserve(m_msh.faces_size());

        const auto num_cell_dofs      = disk::vector_basis_size(m_hdi.cell_degree(), dimension, dimension);
        const auto num_face_dofs      = disk::vector_basis_size(m_hdi.face_degree(), dimension - 1, dimension);
        const auto num_face_dofs_cont = disk::vector_basis_size(m_hdi.face_degree() + 1, dimension - 1, dimension);
        const auto total_dof          = m_msh.cells_size() * num_cell_dofs + m_msh.faces_size() * num_face_dofs;

        for (auto& cl : m_msh)
        {
            const auto ci = disk::contact_info<mesh_type>(m_msh, cl, m_hdi, m_bnd);
            m_solution.push_back(vector_type::Zero(ci.num_total_dofs()));
        }

        for (int i = 0; i < m_msh.faces_size(); i++)
        {
            if (m_bnd.contact_boundary_type(i) == disk::SIGNORINI_FACE)
            {
                m_solution_faces.push_back(vector_type::Zero(num_face_dofs_cont));
            }
            else if (m_bnd.contact_boundary_type(i) == disk::SIGNORINI_CELL)
            {
                m_solution_faces.push_back(vector_type::Zero(0));
            }
            else
            {
                m_solution_faces.push_back(vector_type::Zero(num_face_dofs));
            }
        }

        if (m_verbose)
        {
            std::cout << "Informations about the problem:" << std::endl;
            std::cout << "** Numbers of cells: " << m_msh.cells_size() << std::endl;
            std::cout << "** Numbers of faces: " << m_msh.faces_size()
                      << " ( boundary faces: " << m_msh.boundary_faces_size() << " )" << std::endl;
            std::cout << "** Numbers of dofs: " << total_dof << std::endl;
            std::cout << "** After static condensation: " << std::endl;
            std::cout << "** Numbers of dofs: " << m_msh.faces_size() * num_face_dofs << std::endl;
            std::cout << "** Number of integration points: " << m_law.getNumberOfQP() << std::endl;
            std::cout << " " << std::endl;

            m_law.getMaterialData().print();
        }

        total_dof_depl_static = m_msh.faces_size() * num_face_dofs;
    }

    void
    pre_computation()
    {
        m_gradient_precomputed.clear();
        m_gradient_precomputed.reserve(m_msh.cells_size());

        m_stab_precomputed.clear();
        m_stab_precomputed.reserve(m_msh.cells_size());

        for (auto& cl : m_msh)
        {
            /////// Gradient Reconstruction /////////
            const auto gr = MK::make_matrix_symmetric_gradrec(m_msh, cl, m_hdi, m_bnd);
            m_gradient_precomputed.push_back(gr.first);

            if (m_rp.m_stab)
            {
                switch (m_rp.m_stab_type)
                {
                    case HHO:
                    {
                        const auto recons = MK::make_vector_hho_symmetric_laplacian(m_msh, cl, m_hdi, m_bnd);
                        m_stab_precomputed.push_back(
                          MK::make_vector_hho_stabilization(m_msh, cl, recons.first, m_hdi, m_bnd));
                        break;
                    }
                    case HDG:
                    {
                        m_stab_precomputed.push_back(MK::make_vector_hdg_stabilization(m_msh, cl, m_hdi, m_bnd));
                        break;
                    }
                    // case DG:
                    // {
                    //     m_stab_precomputed.push_back(make_vector_dg_stabilization(m_msh, cl, m_hdi));
                    //     break;
                    // }
                    case NO: { break;
                    }
                    default: throw std::invalid_argument("Unknown stabilization");
                }
            }
        }
    }

    std::list<time_step>
    compute_time_step() const
    {
        // time step
        std::list<time_step> list_step;

        scalar_type prev_time = 0.0;
        for (int n = 0; n < m_rp.m_time_step.size(); n++)
        {
            auto              time_info = m_rp.m_time_step[n];
            const auto        time_curr = time_info.first;
            const scalar_type delta_t   = (time_curr - prev_time) / time_info.second;
            for (int i = 0; i < time_info.second; i++)
            {
                time_step step;
                step.time  = prev_time + (i + 1) * delta_t;
                step.level = 1;
                list_step.push_back(step);
            }
            prev_time = time_curr;
        }

        return list_step;
    }

  public:
    tresca_solver(const mesh_type& msh, const bnd_type& bnd, const param_type& rp, const material_type& material_data) :
      m_msh(msh), m_verbose(rp.m_verbose), m_convergence(false), m_rp(rp), m_bnd(bnd), m_material_data(material_data)
    {
        int face_degree = rp.m_face_degree;
        if (rp.m_face_degree < 0)
        {
            std::cout << "'face_degree' should be > 0. Reverting to 1." << std::endl;
            face_degree = 1;
        }

        m_rp.m_face_degree = face_degree;

        int cell_degree = rp.m_cell_degree;
        if ((face_degree - 1 > cell_degree) or (cell_degree > face_degree + 1))
        {
            std::cout << "'cell_degree' should be 'face_degree + 1' =>"
                      << "'cell_degree' => 'face_degree -1'. Reverting to 'face_degree'." << std::endl;
            cell_degree = face_degree;
        }

        m_rp.m_cell_degree = cell_degree;

        int grad_degree = rp.m_grad_degree;
        if (grad_degree < face_degree)
        {
            std::cout << "'grad_degree' should be >= 'face_degree'. Reverting to 'face_degree'." << std::endl;
            grad_degree = face_degree;
        }

        m_rp.m_grad_degree = grad_degree;

        m_hdi = disk::hho_degree_info(cell_degree, face_degree, grad_degree);

        // compute mesh for post-processing
        gmesh_io = disk::gmsh_io<mesh_type>(m_msh);

        m_law = law_type(m_msh, 2 * m_hdi.grad_degree());
        m_law.addMaterialData(material_data);

        if (m_verbose)
        {
            std::cout << "Informations about the hho's degree:" << std::endl;
            m_hdi.info_degree();
            std::cout << "Informations about the parameters:" << std::endl;
            m_rp.infos();
            std::cout << "Informations about the boundary conditions:" << std::endl;
            m_bnd.boundary_info();
        }

        init();
    }

    bool
    verbose(void) const
    {
        return m_verbose;
    }
    void
    verbose(bool v)
    {
        m_verbose = v;
    }

    template<typename LoadFunction>
    SolverInfo
    compute(const LoadFunction& lf)
    {
        SolverInfo  si;
        timecounter ttot;
        ttot.tic();

        // time step
        std::list<time_step> list_step = compute_time_step();

        int current_step = 0;
        int total_step   = list_step.size();

        scalar_type prev_time = 0.0;

        // time of saving
        bool time_saving(false);
        if (m_rp.m_n_time_save > 0)
        {
            time_saving = true;
        }

        timecounter t1;
        t1.tic();
        this->pre_computation();
        t1.toc();
        if (m_verbose)
            std::cout << "-Precomputation: " << t1.to_double() << " sec" << std::endl;

        // Newton solver
        MK::NewtonRaphson_solver_tresca<mesh_type> newton_solver(m_msh, m_hdi, m_bnd, m_rp);
        newton_solver.initialize(m_solution_faces, m_solution);

        // loading
        while (!list_step.empty())
        {
            current_step += 1;
            const time_step   step         = list_step.front();
            const scalar_type current_time = step.time;

            if (m_verbose)
            {
                std::cout << "------------------------------------------------------------------------"
                             "----------------------"
                          << std::endl;
                std::cout << "************************ Time : " << current_time << " sec (step: " << current_step << "/"
                          << total_step << ", sublevel: " << step.level << " ) *************************|" << std::endl;
            }

            auto rlf = [&lf, &current_time ](const point<scalar_type, dimension>& p) -> auto
            {
                return disk::priv::inner_product(current_time, lf(p));
            };

            m_bnd.multiplyAllFunctionsByAFactor(current_time);

            // correction
            NewtonSolverInfo newton_info =
              newton_solver.compute(rlf, m_gradient_precomputed, m_stab_precomputed, m_law);
            si.updateInfo(newton_info);

            if (m_verbose)
            {
                newton_info.printInfo();
            }

            m_convergence = newton_solver.test_convergence();

            if (!m_convergence)
            {
                if (step.level > m_rp.m_sublevel)
                {
                    std::cout << "***********************************************************" << std::endl;
                    std::cout << "***** PROBLEM OF CONVERGENCE: We stop the calcul here *****" << std::endl;
                    std::cout << "***********************************************************" << std::endl;
                    break;
                }
                else
                {
                    if (m_verbose)
                    {
                        std::cout << "***********************************************************" << std::endl;
                        std::cout << "*****     NO CONVERGENCE: We split the time step     ******" << std::endl;
                        std::cout << "***********************************************************" << std::endl;
                    }
                    total_step += 1;
                    current_step -= 1;
                    time_step new_step;
                    new_step.time  = prev_time + (current_time - prev_time) / 2.0;
                    new_step.level = step.level + 1;
                    list_step.push_front(new_step);
                }
            }
            else
            {
                prev_time = current_time;
                list_step.pop_front();
                newton_solver.save_solutions(m_solution_faces, m_solution);
                m_law.update();

                if (time_saving)
                {
                    if (m_rp.m_time_save.front() < prev_time + 1E-5)
                    {
                        std::cout << "** Save results" << std::endl;
                        std::string name = "result" + std::to_string(dimension) + "D_k" +
                                           std::to_string(m_hdi.face_degree()) + "_l" +
                                           std::to_string(m_hdi.cell_degree()) + "_g" +
                                           std::to_string(m_hdi.grad_degree()) + "_t" + std::to_string(prev_time) + "_";

                        // this->compute_discontinuous_displacement(name + "depl_disc.msh");
                        this->compute_continuous_displacement(name + "depl_cont.msh");
                        // this->compute_discontinuous_stress(name + "stress_disc.msh");
                        // this->compute_continuous_stress(name + "stress_cont.msh");
                        // this->compute_stress_GP(name + "stress_GP.msh");
                        // this->compute_continuous_deformed(name + "deformed_cont.msh");
                        // this->compute_discontinuous_deformed(name + "deformed_disc.msh");

                        m_rp.m_time_save.pop_front();
                        if (m_rp.m_time_save.empty())
                            time_saving = false;
                    }
                }
            }
        }

        si.m_time_step = total_step;

        ttot.toc();
        si.m_time_solver = ttot.to_double();

        return si;
    }

    template<typename As>
    void
    init_solution(const As& as)
    {
        size_t cell_i = 0;

        for (auto& cl : m_msh)
        {
            m_solution.at(cell_i++) = disk::project_function(m_msh, cl, m_hdi, m_bnd, as, 2);
        }

        for (auto itor = m_msh.faces_begin(); itor != m_msh.faces_end(); itor++)
        {
            const auto bfc     = *itor;
            const auto face_id = m_msh.lookup(bfc);

            if (m_bnd.contact_boundary_type(face_id) == disk::SIGNORINI_FACE)
            {
                const auto proj_bcf = disk::project_function(m_msh, bfc, m_hdi.face_degree() + 1, as, 2);
                assert(m_solution_faces[face_id].size() == proj_bcf.size());

                m_solution_faces[face_id] = proj_bcf;
            }
            else if (m_bnd.contact_boundary_type(face_id) == disk::SIGNORINI_CELL)
            {
                assert(m_solution_faces[face_id].size() == 0);
            }
            else
            {
                const auto proj_bcf = disk::project_function(m_msh, bfc, m_hdi.face_degree(), as, 2);
                assert(m_solution_faces[face_id].size() == proj_bcf.size());

                m_solution_faces[face_id] = proj_bcf;
            }
        }
    }

    bool
    test_convergence() const
    {
        return m_convergence;
    }

    int
    getDofs()
    {
        return total_dof_depl_static;
    }

    void
    printSolutionCell() const
    {
        int        cell_i = 0;
        const auto cbs    = disk::vector_basis_size(m_hdi.cell_degree(), dimension, dimension);

        std::cout << "Solution at the cells:" << std::endl;
        for (auto& cl : m_msh)
        {
            std::cout << "cell " << cell_i << ": " << std::endl;
            std::cout << m_solution.at(cell_i++).head(cbs) << std::endl;
        }
    }

    // compute l2 error
    template<typename AnalyticalSolution>
    scalar_type
    compute_l2_displacement_error(const AnalyticalSolution& as)
    {
        scalar_type err_dof = 0;

        const auto cbs = disk::vector_basis_size(m_hdi.cell_degree(), dimension, dimension);

        int cell_i = 0;

        for (auto& cl : m_msh)
        {
            const auto x = m_solution.at(cell_i++).head(cbs);

            const vector_type true_dof = disk::project_function(m_msh, cl, m_hdi.cell_degree(), as, 2);

            const auto        cb   = disk::make_vector_monomial_basis(m_msh, cl, m_hdi.cell_degree());
            const matrix_type mass = disk::make_mass_matrix(m_msh, cl, cb);

            const vector_type comp_dof = x.head(cbs);
            // std::cout << "cell:" << std::endl;
            // std::cout << true_dof.transpose() << std::endl;
            // std::cout << comp_dof.transpose() << std::endl;
            const vector_type diff_dof = (true_dof - comp_dof);
            assert(comp_dof.size() == true_dof.size());
            err_dof += diff_dof.dot(mass * diff_dof);
        }

        return sqrt(err_dof);
    }

    // compute l2 error
    template<typename AnalyticalSolution>
    scalar_type
    compute_H1_error(const AnalyticalSolution& as)
    {
        scalar_type err_dof = 0;

        int cell_i = 0;

        for (auto& cl : m_msh)
        {
            const auto gr   = MK::make_matrix_symmetric_gradrec(m_msh, cl, m_hdi, m_bnd);
            const auto stab = m_stab_precomputed.at(cell_i);
            // const auto recons = make_vector_hho_symmetric_laplacian(m_msh, cl, m_hdi);
            // const auto stab   = make_vector_hho_stabilization(m_msh, cl, recons.first, m_hdi);
            // const auto stab = MK::make_vector_hdg_stabilization(m_msh, cl, m_hdi, m_bnd);
            const auto Ah = gr.second + stab;

            const auto x = m_solution.at(cell_i++);

            const vector_type true_dof = disk::project_function(m_msh, cl, m_hdi, m_bnd, as, 2);

            // std::cout << "cell:" << std::endl;
            // std::cout << true_dof.transpose() << std::endl;
            // std::cout << x.transpose() << std::endl;
            const vector_type diff_dof = (true_dof - x);

            err_dof += diff_dof.dot(Ah * diff_dof);
        }

        return sqrt(err_dof);
    }

    // compute l2 error
    template<typename AnalyticalSolution>
    scalar_type
    compute_l2_pres_error(const AnalyticalSolution& as)
    {
        scalar_type err_dof = 0;

        int cell_i = 0;

        MK::tresca<mesh_type> elem(m_msh, m_hdi, m_material_data, m_rp, m_bnd);

        for (auto& cl : m_msh)
        {
            const auto gb           = disk::make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());
            const auto law_cell     = m_law.getCellIVs(cell_i);
            const auto stress_coeff = law_cell.projectStressSymOnCell(m_msh, cl, m_hdi, m_material_data);

            const auto hT = diameter(m_msh, cl);

            const auto fcs    = m_bnd.faces_with_contact(cl);
            size_t     face_i = 0;
            for (auto& fc : fcs)
            {
                const auto fb           = disk::make_vector_monomial_basis(m_msh, fc, m_hdi.face_degree());
                const auto n            = normal(m_msh, cl, fc);
                const auto qp_deg       = std::max(m_hdi.cell_degree(), m_hdi.grad_degree());
                const auto qps          = integrate(m_msh, fc, 2 * qp_deg + 4);
                const auto hF           = diameter(m_msh, fc);
                const auto gamma_F      = m_rp.m_gamma_0 / hF;
                const auto contact_type = m_bnd.contact_boundary_type(fc);
                const auto contact_id   = m_bnd.contact_boundary_id(fc);

                // if(contact_id == 1)
                // {
                for (auto& qp : qps)
                {
                    const scalar_type pres_sol = as(qp.point());

                    const auto stress = elem.eval_stress(stress_coeff, gb, qp.point());

                    const scalar_type pres_comp = stress(1, 1);

                    // std::cout << pres_sol << " vs " << pres_comp << std::endl;

                    err_dof += hT * qp.weight() * (pres_sol - pres_comp) * (pres_sol - pres_comp);
                }
                // }
            }
        }

        return sqrt(err_dof);
    }

    void
    compute_discontinuous_displacement(const std::string& filename) const
    {
        gmesh_io.save_vector_hho_solution_discontinuous(filename, m_msh, m_hdi, m_solution);
    }

    void
    compute_continuous_displacement(const std::string& filename) const
    {
        gmesh_io.save_vector_hho_solution_continuous(filename, m_msh, m_hdi, m_solution);
    }

    void
    compute_continuous_deformed(const std::string& filename) const
    {
        gmesh_io.save_vector_hho_deformed_continuous(filename, m_msh, m_hdi, m_solution);
    }

    void
    compute_discontinuous_deformed(const std::string& filename) const
    {
        gmesh_io.save_vector_hho_deformed_discontinuous(filename, m_msh, m_hdi, m_solution);
    }

    void
    compute_stress_GP(const std::string& filename) const
    {
        const gmsh::Gmesh& gmsh = gmesh_io.gmesh();

        std::vector<gmsh::Data>    data;    // create data (not used)
        std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point
        size_t                     nb_nodes(gmsh.getNumberofNodes());

        const auto material_data = m_law.getMaterialData();

        int cell_i = 0;
        for (auto& cl : m_msh)
        {
            const auto law_quadpoints = m_law.getCellIVs(cell_i).getIVs();

            // Loop on nodes
            for (auto& qp : law_quadpoints)
            {
                const auto                stress = qp.compute_stress(material_data);
                const std::vector<double> tens   = disk::convertToVectorGmsh(stress);

                // Add GP
                // Create a node at gauss point
                nb_nodes++;
                const gmsh::Node    new_node = disk::convertPoint(qp.point(), nb_nodes);
                const gmsh::SubData sdata(tens, new_node);
                subdata.push_back(sdata); // add subdata
            }
            cell_i++;
        }

        // Save
        gmsh::NodeData nodedata(9, 0.0, "Stress_GP", data, subdata); // create and init a nodedata view

        nodedata.saveNodeData(filename, gmsh); // save the view
    }

    void
    compute_discontinuous_stress(const std::string& filename) const
    {
        gmsh::Gmesh gmsh(dimension);
        auto        storage = m_msh.backend_storage();

        std::vector<gmsh::Data>          data;    // create data (not used)
        const std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point

        const size_t                grad_degree = m_hdi.grad_degree();
        const MK::tresca<mesh_type> elem(m_msh, m_hdi, m_material_data, m_rp, m_bnd);

        int    cell_i   = 0;
        size_t nb_nodes = 0;
        for (auto& cl : m_msh)
        {
            const auto              gb           = disk::make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());
            const auto              law_cell     = m_law.getCellIVs(cell_i);
            const auto              stress_coeff = law_cell.projectStressSymOnCell(m_msh, cl, m_hdi, m_material_data);
            const auto              cell_nodes   = disk::cell_nodes(m_msh, cl);
            std::vector<gmsh::Node> new_nodes;

            // loop on the nodes of the cell
            for (int i = 0; i < cell_nodes.size(); i++)
            {
                nb_nodes++;
                const auto point_ids = cell_nodes[i];
                const auto pt        = storage->points[point_ids];

                const auto stress = elem.eval_stress(stress_coeff, gb, pt);

                const std::vector<double>   tens = disk::convertToVectorGmsh(stress);
                const std::array<double, 3> coor = disk::init_coor(pt);

                // Add a node
                const gmsh::Node tmp_node(coor, nb_nodes, 0);
                new_nodes.push_back(tmp_node);
                gmsh.addNode(tmp_node);

                const gmsh::Data datatmp(nb_nodes, tens);
                data.push_back(datatmp);
            }

            // Add new element
            disk::add_element(gmsh, new_nodes);
            cell_i++;
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(9, 0.0, "stress_node_disc", data, subdata);

        // Save the view
        nodedata.saveNodeData(filename, gmsh);
    }

    void
    compute_continuous_stress(const std::string& filename) const
    {
        const gmsh::Gmesh& gmsh      = gmesh_io.gmesh();
        const auto&        post_mesh = gmesh_io.post_mesh();
        auto               storage   = post_mesh.mesh().backend_storage();

        const static_matrix<scalar_type, dimension, dimension> vzero =
          static_matrix<scalar_type, dimension, dimension>::Zero();

        const size_t nb_nodes(gmsh.getNumberofNodes());

        // first(number of data at this node), second(cumulated value)
        std::vector<std::pair<size_t, static_matrix<scalar_type, dimension, dimension>>> value(
          nb_nodes, std::make_pair(0, vzero));

        const size_t                grad_degree = m_hdi.grad_degree();
        const MK::tresca<mesh_type> elem(m_msh, m_hdi, m_material_data, m_rp, m_bnd);

        int cell_i = 0;

        for (auto& cl : m_msh)
        {
            const auto gb           = disk::make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());
            const auto law_cell     = m_law.getCellIVs(cell_i);
            const auto stress_coeff = law_cell.projectStressSymOnCell(m_msh, cl, m_hdi, m_material_data);
            const auto cell_nodes   = post_mesh.nodes_cell(cell_i);

            // Loop on the nodes of the cell
            for (auto& point_id : cell_nodes)
            {
                const auto pt = storage->points[point_id];

                const auto stress = elem.eval_stress(stress_coeff, gb, pt);

                value[point_id].first++;
                value[point_id].second += stress;
            }
            cell_i++;
        }

        std::vector<gmsh::Data>    data;    // create data
        std::vector<gmsh::SubData> subdata; // create subdata
        data.reserve(nb_nodes);             // data has a size of nb_node

        // Compute the average value and save it
        for (int i_node = 0; i_node < value.size(); i_node++)
        {
            const static_matrix<scalar_type, dimension, dimension> stress_avr =
              value[i_node].second / double(value[i_node].first);

            const gmsh::Data tmp_data(i_node + 1, disk::convertToVectorGmsh(stress_avr));
            data.push_back(tmp_data);
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(9, 0.0, "stress_node_cont", data, subdata);
        // Save the view
        nodedata.saveNodeData(filename, gmsh);
    }

    void
    compute_equivalent_plastic_strain_GP(const std::string& filename) const
    {
        const gmsh::Gmesh& gmsh = gmesh_io.gmesh();

        std::vector<gmsh::Data>    data;    // create data (not used)
        std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point
        size_t                     nb_nodes(gmsh.getNumberofNodes());

        const auto material_data = m_law.getMaterialData();

        int cell_i = 0;
        for (auto& cl : m_msh)
        {

            const auto law_quadpoints = m_law.getCellIVs(cell_i).getIVs();

            // Loop on nodes
            for (auto& qp : law_quadpoints)
            {
                const auto                p   = qp.getAccumulatedPlasticStrain();
                const std::vector<double> p_s = disk::convertToVectorGmsh(p);

                // Add GP
                // Create a node at gauss point
                nb_nodes++;
                const gmsh::Node    new_node = disk::convertPoint(qp.point(), nb_nodes);
                const gmsh::SubData sdata(p_s, new_node);
                subdata.push_back(sdata); // add subdata
            }
            cell_i++;
        }

        // Save
        gmsh::NodeData nodedata(1, 0.0, "p_GP", data, subdata); // create and init a nodedata view

        nodedata.saveNodeData(filename, gmsh); // save the view
    }

    void
    compute_is_plastic_GP(const std::string& filename) const
    {
        gmsh::Gmesh gmsh = gmesh_io.gmesh();

        std::vector<gmsh::Data>    data;    // create data (not used)
        std::vector<gmsh::SubData> subdata; // create subdata to save soution at gauss point
        size_t                     nb_nodes(gmsh.getNumberofNodes());

        const auto material_data = m_law.getMaterialData();

        int cell_i = 0;
        for (auto& cl : m_msh)
        {

            const auto law_quadpoints = m_law.getCellIVs(cell_i).getIVs();

            // Loop on nodes
            for (auto& qp : law_quadpoints)
            {
                scalar_type p = 0;
                if (qp.is_plastic())
                    p = 1;

                const std::vector<double> p_s = disk::convertToVectorGmsh(p);

                // Add GP
                // Create a node at gauss point
                nb_nodes++;
                const gmsh::Node    new_node = disk::convertPoint(qp.point(), nb_nodes);
                const gmsh::SubData sdata(p_s, new_node);
                subdata.push_back(sdata); // add subdata
            }
            cell_i++;
        }

        // Save
        gmsh::NodeData nodedata(1, 0.0, "state_GP", data, subdata); // create and init a nodedata view

        nodedata.saveNodeData(filename, gmsh); // save the view
    }

    void
    compute_continuous_pres(const std::string& filename) const
    {
        const gmsh::Gmesh& gmsh      = gmesh_io.gmesh();
        const auto&        post_mesh = gmesh_io.post_mesh();
        auto               storage   = post_mesh.mesh().backend_storage();

        const static_matrix<scalar_type, dimension, dimension> vzero =
          static_matrix<scalar_type, dimension, dimension>::Zero();

        const size_t nb_nodes(gmsh.getNumberofNodes());

        // first(number of data at this node), second(cumulated value)
        std::vector<std::pair<size_t, scalar_type>> value(nb_nodes, std::make_pair(0, 0.0));

        const size_t                grad_degree = m_hdi.grad_degree();
        const MK::tresca<mesh_type> elem(m_msh, m_hdi, m_material_data, m_rp, m_bnd);

        int cell_i = 0;

        for (auto& cl : m_msh)
        {
            const auto        fcs = disk::faces(m_msh, cl);
            std::vector<bool> fcs_bnd;
            bool              at_least_one = false;

            for (auto fc : fcs)
            {
                const bool l_bnd = m_msh.is_boundary(fc);
                fcs_bnd.push_back(l_bnd);

                if (l_bnd)
                {
                    at_least_one = true;
                }
            }

            if (at_least_one)
            {
                const auto gb           = disk::make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());
                const auto law_cell     = m_law.getCellIVs(cell_i);
                const auto stress_coeff = law_cell.projectStressSymOnCell(m_msh, cl, m_hdi, m_material_data);

                size_t face_i = 0;
                for (auto fc : fcs)
                {

                    if (fcs_bnd[face_i++])
                    {
                        const auto face_nodes = post_mesh.nodes_face(m_msh.lookup(fc));
                        const auto n          = disk::normal(m_msh, cl, fc);

                        // Loop on the nodes of the cell
                        for (auto& point_id : face_nodes)
                        {
                            const auto pt = storage->points[point_id];

                            const auto stress_nn = elem.eval_stress_nn(stress_coeff, gb, n, pt);

                            value[point_id].first++;
                            value[point_id].second += stress_nn;
                        }
                    }
                }
            }

            cell_i++;
        }

        std::vector<gmsh::Data>    data;    // create data
        std::vector<gmsh::SubData> subdata; // create subdata
        data.reserve(nb_nodes);             // data has a size of nb_node

        // Compute the average value and save it
        for (int i_node = 0; i_node < value.size(); i_node++)
        {
            scalar_type stress_avr = 0.0;

            if (value[i_node].first > 0)
                stress_avr = value[i_node].second / double(value[i_node].first);

            const gmsh::Data tmp_data(i_node + 1, disk::convertToVectorGmsh(stress_avr));
            data.push_back(tmp_data);
        }

        // Create and init a nodedata view
        gmsh::NodeData nodedata(1, 0.0, "press_node_cont", data, subdata);
        // Save the view
        nodedata.saveNodeData(filename, gmsh);
    }

    void
    contact_quantities(const std::string& filename) const
    {
        std::ofstream output;
        output.open(filename, std::ofstream::out | std::ofstream::trunc);

        if (!output.is_open())
        {
            std::cerr << "Unable to open file " << filename << std::endl;
        }

        output << "#x"
               << "\t"
               << "y"
               << "\t"
               << "u_n"
               << "\t"
               << "sigma_nn"
               << "\t"
               << "phi_n_1"
               << "\t"
               << "[phi_n_1]-"
               << "\t"
               << "u_T"
               << "\t"
               << "sigma_nt"
               << "\t"
               << "phi_t_1"
               << "\t"
               << "[phi_n_1]s"
               << "\t"
               << "abs(phi_t_1)"
               << "\t"
               << "abs([phi_n_1]s)"
               << "\t"
               << "abs([phi_n_1]s)/s"
               << "\t"
               << "r" << std::endl;

        const size_t grad_degree = m_hdi.grad_degree();
        const auto   cbs         = disk::vector_basis_size(m_hdi.cell_degree(), dimension, dimension);

        MK::tresca<mesh_type> elem(m_msh, m_hdi, m_material_data, m_rp, m_bnd);

        int cell_i = 0;
        for (auto& cl : m_msh)
        {
            if (m_bnd.cell_has_contact_faces(cl))
            {
                const auto ci = disk::contact_info<mesh_type>(m_msh, cl, m_hdi, m_bnd);

                const auto gb           = disk::make_sym_matrix_monomial_basis(m_msh, cl, m_hdi.grad_degree());
                const auto law_cell     = m_law.getCellIVs(cell_i);
                const auto stress_coeff = law_cell.projectStressSymOnCell(m_msh, cl, m_hdi, m_material_data);

                const auto uTF = m_solution.at(cell_i);

                const auto cb = disk::make_vector_monomial_basis(m_msh, cl, m_hdi.cell_degree());

                const auto fcs = ci.faces();

                size_t offset = cbs;
                for (auto& fc : fcs)
                {
                    const auto facedeg = ci.face_degree(m_bnd, fc);
                    const auto fbs     = ci.num_face_dofs(m_bnd, fc);

                    if (m_bnd.is_contact_face(fc))
                    {
                        const auto fb     = disk::make_vector_monomial_basis(m_msh, fc, facedeg);
                        const auto n      = normal(m_msh, cl, fc);
                        const auto qp_deg = 0;
                        // std::max(m_hdi.cell_degree(), m_hdi.face_degree());
                        const auto qps          = integrate(m_msh, fc, 2 * qp_deg);
                        const auto hF           = diameter(m_msh, fc);
                        const auto gamma_F      = m_rp.m_gamma_0 / hF;
                        const auto contact_type = m_bnd.contact_boundary_type(fc);

                        const auto s_func = m_bnd.contact_boundary_func(fc);

                        for (auto& qp : qps)
                        {
                            const scalar_type sigma_nn_u = elem.eval_stress_nn(stress_coeff, gb, n, qp.point());
                            const auto        sigma_nt_u = elem.eval_stress_nt(stress_coeff, gb, n, qp.point());

                            const scalar_type r =
                              std::sqrt(qp.point().x() * qp.point().x() + qp.point().y() * qp.point().y());

                            if (contact_type == disk::SIGNORINI_CELL)
                            {
                                const scalar_type phi_n_1_u =
                                  elem.eval_phi_n_uT(stress_coeff, gb, cb, uTF, n, gamma_F, qp.point());
                                const scalar_type phi_n_1_u_proj =
                                  elem.eval_proj_phi_n_uT(stress_coeff, gb, cb, uTF, n, gamma_F, qp.point());
                                const scalar_type uT_n_u = elem.eval_uT_n(cb, uTF, n, qp.point());
                                const auto        phi_t_1_u =
                                  elem.eval_phi_t_uT(stress_coeff, gb, cb, uTF, n, gamma_F, qp.point());
                                const auto phi_t_1_u_proj = elem.eval_proj_phi_t_uT(
                                  stress_coeff, gb, cb, uTF, n, gamma_F, s_func(qp.point()), qp.point());

                                const auto uT_t_u = elem.eval_uT_t(cb, uTF, n, qp.point());

                                output << qp.point().x() << "\t" << qp.point().y() << "\t" << uT_n_u << "\t"
                                       << sigma_nn_u << "\t" << phi_n_1_u << "\t" << phi_n_1_u_proj << "\t"
                                       << uT_t_u.transpose() << "\t" << sigma_nt_u.transpose() << "\t"
                                       << phi_t_1_u.transpose() << "\t" << phi_t_1_u_proj.transpose() << "\t"
                                       << phi_t_1_u.norm() << "\t" << phi_t_1_u_proj.norm() << "\t"
                                       << phi_t_1_u_proj.norm() / m_rp.m_threshold << "\t" << r << std ::endl;
                            }
                            else
                            {
                                const auto        uF = uTF.segment(offset, fbs);
                                const scalar_type phi_n_1_u =
                                  elem.eval_phi_n_uF(fc, stress_coeff, gb, fb, uTF, offset, n, gamma_F, qp.point());
                                const scalar_type phi_n_1_u_proj = elem.eval_proj_phi_n_uF(
                                  fc, stress_coeff, gb, fb, uTF, offset, n, gamma_F, qp.point());
                                const scalar_type uT_n_u = elem.eval_uF_n(fb, uF, n, qp.point());
                                const auto        phi_t_1_u =
                                  elem.eval_phi_t_uF(stress_coeff, gb, fb, uTF, offset, n, gamma_F, qp.point());
                                const auto phi_t_1_u_proj = elem.eval_proj_phi_t_uF(
                                  stress_coeff, gb, fb, uTF, offset, n, gamma_F, s_func(qp.point()), qp.point());
                                const auto uF_t_u = elem.eval_uF_t(fb, uF, n, qp.point());

                                output << qp.point().x() << "\t" << qp.point().y() << "\t" << uT_n_u << "\t"
                                       << sigma_nn_u << "\t" << phi_n_1_u << "\t" << phi_n_1_u_proj << "\t"
                                       << uF_t_u.transpose() << "\t" << sigma_nt_u.transpose() << "\t"
                                       << phi_t_1_u.transpose() << "\t" << phi_t_1_u_proj.transpose() << "\t"
                                       << phi_t_1_u.norm() << "\t" << phi_t_1_u_proj.norm() << "\t"
                                       << phi_t_1_u_proj.norm() / m_rp.m_threshold << "\t" << r << std ::endl;
                            }
                        }
                    }

                    offset += fbs;
                }
            }
            cell_i++;
        }
        output.close();
    }
};