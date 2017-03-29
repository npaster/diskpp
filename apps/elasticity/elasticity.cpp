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

#include <iostream>

#include "loaders/loader.hpp"
#include "elasticity.hpp"
#include "timecounter.h"
#include "assert.h"

#define _USE_MATH_DEFINES
#include <cmath>

template<typename MeshType, typename Function, typename Solution>
void test_new_elasticity(MeshType& msh, const Function& load, const Solution& solution, size_t degree)
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
                                       cell_basis_type, cell_quadrature_type,
                                       face_basis_type, face_quadrature_type> gradrec(degree);



   disk::elas_like_stabilization<mesh_type,
                                 cell_basis_type, cell_quadrature_type,
                                 face_basis_type, face_quadrature_type> stab(degree);

   disk::diffusion_like_static_condensation<mesh_type,
                                             cell_basis_type,
                                             cell_quadrature_type,
                                             face_basis_type,
                                             face_quadrature_type> statcond(degree);

   disk::assembler_elas<mesh_type, face_basis_type, face_quadrature_type> assembler(msh, degree);


   disk::l2_error<mesh_type, cell_basis_type, cell_quadrature_type, face_basis_type,
   face_quadrature_type, Solution> error(degree);

   disk::projector_elas<mesh_type,
                           cell_basis_type, cell_quadrature_type,
                           face_basis_type, face_quadrature_type> proj(degree);

   scalar_type mu      = 1.0;
   scalar_type lambda  = 0.0;
   const size_t DIM = msh.dimension;

   timecounter tc;

   tc.tic();
   size_t cell_i = 0;
   for (auto& cl : msh)
   {
      gradrec.compute(msh, cl);
      //          std::cout << "gradrec" << std::endl;
      //          std::cout << gradrec.oper << std::endl;




      stab.compute(msh, cl, gradrec.oper);

      //          std::cout << "stab" << std::endl;
      //          std::cout << stab.data << std::endl;

      //std::cout << "test gradrec" << std::endl;
//       test_gradient_reconstruction<mesh_type, cell_basis_type, cell_quadrature_type,
//       face_basis_type, face_quadrature_type, Solution>
//       (msh, cl, solution, gradrec.oper, degree);

      auto dof_true = proj.compute_whole(msh, cl, solution);

      //   std::cout << "STU true" << std::endl;
      //   std::cout << stab.data * dof_true << std::endl;
      // std::cout << "USTU true " << dof_true.transpose() * stab.data * dof_true << std::endl;

      auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, load, degree);


      /*
       s td::cout << "rhs" << std::endl;            *
       std::cout << cell_rhs << std::endl;*/

      assert(cell_rhs.size() == DIM * binomial(degree + DIM, degree));
      dynamic_matrix<scalar_type> loc =   2.0* mu * gradrec.data +
                                          2.0*mu * stab.data;
      auto sc = statcond.compute(msh, cl, loc, cell_rhs);

      //         std::cout << "mat_cond" << std::endl;
      //         std::cout << sc.first << std::endl;
      //         std::cout << "rhs_cond" << std::endl;
      //         std::cout << sc.second << std::endl;


      assembler.assemble(msh, cl, sc);
      cell_i++;
      //          if(cell_i == 4)
      //             assert(false);
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
    s td::cout << "* * Solving for " << systsz << " unknow*ns." << std::endl;
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
      dynamic_matrix<scalar_type> loc =   2.0* mu * gradrec.data +
                                          2.0* mu * stab.data;

      dynamic_vector<scalar_type> x = statcond.recover(msh, cl, loc, cell_rhs, xFs);

      auto dof_true = proj.compute_whole(msh, cl, solution);

      //   std::cout << "STU" << std::endl;
      //   std::cout << stab.data * x << std::endl;
      //std::cout << "USTU " << x.transpose() * stab.data * x << std::endl;

      //   std::cout << "STU true" << std::endl;
      //   std::cout << stab.data * dof_true << std::endl;
      //std::cout << "USTU true " << dof_true.transpose() * stab.data * dof_true << std::endl;

      //         std::cout << "recover" << std::endl;
      //         std::cout << x << std::endl;

      l2 += error.compute(msh, cl, x, solution);
      cont_cell += 1;

   }

   std::cout << "elasticity l2_error: " << l2 << std::endl;
}



int main(void)
{
   using RealType = double;

   char    *mesh_filename  = nullptr;


   typedef disk::simplicial_mesh<RealType, 2>           mesh_type;




   mesh_type msh;
   disk::netgen_mesh_loader<RealType, 2> loader;
   //if (!loader.read_mesh("/home/C00976/Documents/Disk++/meshes/3D_tetras/netgen/convt01.mesh"))
   if (!loader.read_mesh("/home/C00976/Documents/Disk++/meshes/2D_triangles/netgen/tri01.mesh2d"))
   //if (!loader.read_mesh("/users/npignet/Documents/Diskpp/meshes/2D_triangles/netgen/tri01.mesh2d"))
   {
      std::cout << "Problem loading mesh." << std::endl;
      return 1;
   }
   loader.populate_mesh(msh);

   //auto msh = disk::load_netgen_2d_mesh<RealType>("/home/C00976/Documents/Disk++/meshes/2D_triangles/netgen/tri01.mesh2d");

   typedef RealType scalar_type;

   typedef static_vector<scalar_type, 2> result_type;



   auto f2lin = [](const point<scalar_type,2>& p) -> result_type {
      const scalar_type lambda =1.0;
      const scalar_type mu = 1.0;
      scalar_type fx = 2.*M_PI*M_PI*sin(M_PI*p.x())*sin(M_PI*p.y());

      scalar_type fy = 2.*M_PI*M_PI*cos(M_PI*p.x())*cos(M_PI*p.y());


      return result_type{0.0,0.0};
   };

   auto sf2lin = [](const point<scalar_type,2>& p) -> result_type {
      const scalar_type lambda =1.0;
      const scalar_type mu =  1.0;
      scalar_type fx = p.x();
      scalar_type fy = p.y();

      return result_type{fx,fy};
   };

   auto f2 = [](const point<scalar_type,2>& p) -> result_type {
      const scalar_type lambda =1.0;
      const scalar_type mu = 1.0;
      scalar_type fx = 2.*M_PI*M_PI*sin(M_PI*p.x())*sin(M_PI*p.y());

      scalar_type fy = 2.*M_PI*M_PI*cos(M_PI*p.x())*cos(M_PI*p.y());


      return result_type{fx,fy};
   };

   auto sf2 = [](const point<scalar_type,2>& p) -> result_type {
      const scalar_type lambda =1.0;
      const scalar_type mu =  1.0;
      scalar_type fx = sin(M_PI*p.x())*sin(M_PI*p.y()) /*+ 1.0/(2.0*lambda) * p.x() */;
      scalar_type fy = cos(M_PI*p.x())*cos(M_PI*p.y()) /*+ 1.0/(2.0*lambda) * p.y() */;

      return result_type{fx,fy};
   };






   //     auto f = [](const point<scalar_type,3>& p) -> result_type {
   //         const scalar_type lambda =1.0;
   //         scalar_type fx = M_PI*M_PI*(12*lambda*cos(2*M_PI*p.x())*sin(2*M_PI*p.y())*sin(2*M_PI*p.z())
   //                + 8 * (3*cos(2*M_PI*p.x())-1)*sin(2*M_PI*p.y())*sin(2*M_PI*p.z())
   //                - cos(M_PI*p.x())*sin(M_PI*(p.y()+p.z()))
   //                + (1 + 3./(1+lambda))* sin(M_PI*p.x())*sin(M_PI*p.y())*sin(M_PI*p.z()) );
   //
   //         scalar_type fy = M_PI*M_PI*(12*lambda*cos(2*M_PI*p.y())*sin(2*M_PI*p.x())*sin(2*M_PI*p.z())
   //                + 8 * (3*cos(2*M_PI*p.y())-1)*sin(2*M_PI*p.x())*sin(2*M_PI*p.z())
   //                - cos(M_PI*p.y())*sin(M_PI*(p.x()+p.z())) +
   //                (1 + 3./(1+lambda))* sin(M_PI*p.x())*sin(M_PI*p.y())*sin(M_PI*p.z()) );
   //
   //         scalar_type fz = M_PI*M_PI*(12*lambda*cos(2*M_PI*p.z())*sin(2*M_PI*p.y())*sin(2*M_PI*p.x())
   //                + 8 * (3*cos(2*M_PI*p.z())-1)*sin(2*M_PI*p.y())*sin(2*M_PI*p.x())
   //                - cos(M_PI*p.z())*sin(M_PI*(p.y()+p.x())) +
   //                (1 + 3./(1+lambda))* sin(M_PI*p.x())*sin(M_PI*p.y())*sin(M_PI*p.z()) );
   //
   //         return result_type{fx,fy,fz};
   //     };
   //
   //     auto sf = [](const point<scalar_type,3>& p) -> result_type {
   //         const scalar_type lambda =1.0;
   //         scalar_type fx = sin(2*M_PI*p.y())*sin(2*M_PI*p.z())*(-1 + cos(2*M_PI*p.x()))
   //                         + (1./(1+lambda))* sin(M_PI*p.x())*sin(M_PI*p.y())*sin(M_PI*p.z());
   //         scalar_type fy = sin(2*M_PI*p.z())*sin(2*M_PI*p.x())*(-1 + cos(2*M_PI*p.y()))
   //                         + (1./(1+lambda))* sin(M_PI*p.x())*sin(M_PI*p.y())*sin(M_PI*p.z());
   //         scalar_type fz = sin(2*M_PI*p.x())*sin(2*M_PI*p.y())*(-1 + cos(2*M_PI*p.z()))
   //                         + (1./(1+lambda))* sin(M_PI*p.x())*sin(M_PI*p.y())*sin(M_PI*p.z());
   //         return result_type{fx,fy,fz};
   //     };


   //test_condensation(msh, sf2, 1);
//    test_new_elasticity(msh, f2, sf2, 1);
//    test_new_elasticity(msh, f2, sf2, 2);
//    test_new_elasticity(msh, f2, sf2, 4);
   //test_new_elasticity(msh, f2, sf2, 6);




   ////////////////// 1D //////////////////////////
   std::cout << "Mesh format: 1D uniform" << std::endl;
   int elems_1d = 2;
   auto msh1D = disk::load_uniform_1d_mesh<RealType>(0, 1, elems_1d);
   //test_condensation(msh1D, sf1, 1);

   auto load = [](const point<scalar_type,1>& p) -> auto {
      const scalar_type lambda = scalar_type{0.0};
      const scalar_type mu = scalar_type{1.0};
      return 0.0;
   };

   auto solution = [](const point<scalar_type,1>& p) -> auto {
      return p.x();
   };

   auto f1 = [](const point<scalar_type, 1>& p) -> auto {
      return  2.0*M_PI *  M_PI * sin(p.x() * M_PI);
   };

   auto sf1 = [](const point<scalar_type, 1>& p) -> auto {
      return sin(p.x() * M_PI);
   };

   std::cout << "test 1D " << std::endl;
//    test_gradient(msh1D, f1, sf1, 2);
//    test_gradient(msh1D, f1, sf1, 4);
//    test_gradient(msh1D, f1, sf1, 6);


/*
   auto fe1 = [](const point<scalar_type, 1>& p) -> auto {
      return  3.0*M_PI *  M_PI * sin(p.x() * M_PI);
   };

   auto sfe1 = [](const point<scalar_type, 1>& p) -> auto {
      return sin(p.x() * M_PI);
   };
   //    test_new_elasticity(msh1D, fe1, sfe1, 1);
   //      test_new_elasticity(msh1D, fe1, sfe1, 3);
   //      test_new_elasticity(msh1D, fe1, sfe1, 4);
   //      test_new_elasticity(msh1D, fe1, sfe1, 5);*/


   std::cout << "test 2D " << std::endl;

   auto sfgrad2 = [](const point<scalar_type,2>& p) -> result_type {
      const scalar_type mu =  1.0;
      scalar_type fx = sin(M_PI*p.x()) * sin(M_PI*p.y()) ;
      scalar_type fy = cos(M_PI*p.x()) * cos(M_PI*p.y());

      return result_type{fx,fy};
   };

   auto fgrad2 = [](const point<scalar_type,2>& p) -> result_type {
      const scalar_type mu = 1.0;
      scalar_type fx = 4.*mu*M_PI*M_PI*sin(M_PI*p.x())*sin(M_PI*p.y());
      scalar_type fy = 4.*mu*M_PI*M_PI*cos(M_PI*p.x())*cos(M_PI*p.y());

      return result_type{fx,fy};
   };


   test_gradient(msh, fgrad2, sfgrad2, 1);
   test_gradient(msh, fgrad2, sfgrad2, 3);
   test_gradient(msh, fgrad2, sfgrad2, 6);

   return 0;

}
