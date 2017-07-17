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
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "colormanip.h"

#include "../../config.h"

#ifdef HAVE_SOLVER_WRAPPERS
#include "agmg/agmg.hpp"
#endif

#include "loaders/loader.hpp"
#include "hho/hho.hpp"
#include "geometry/geometry.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "common/eigen.hpp"
#include "hho/hho_bq.hpp"
#include "hho/hho_vector.hpp"



template< typename mesh_type, typename T>
void
test_grad_vector(const mesh_type& msh, const size_t degree, const T epsilon)
{
   typedef typename mesh_type::scalar_type            scalar_type;
   typedef typename mesh_type::cell                   cell_type;
   typedef typename mesh_type::face                   face_type;

   typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
   typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;
   typedef disk::quadrature<mesh_type, cell_type>      matrix_quadrature_type;

   typedef disk::scaled_monomial_vector_basis<mesh_type, cell_type>    cell_basis_type;
   typedef disk::scaled_monomial_vector_basis<mesh_type, face_type>    face_basis_type;
   typedef disk::scaled_monomial_matrix_basis<mesh_type, cell_type>    matrix_basis_type;

   typedef disk::basis_quadrature_data_full<mesh_type, disk::scaled_monomial_vector_basis,
   disk::scaled_monomial_matrix_basis,
   disk::quadrature> bqdata_type;

   typedef disk::gradient_reconstruction_elas_full_bq<bqdata_type>     gradrecopt_type;


   typedef disk::gradient_reconstruction_elas_full< mesh_type,
   cell_basis_type, cell_quadrature_type,
   face_basis_type, face_quadrature_type,
   matrix_basis_type, matrix_quadrature_type>          gradrec_type;



   typedef dynamic_matrix<scalar_type>         matrix_dynamic;

   bqdata_type bqd(degree, degree, degree);

   gradrecopt_type gradrec_opt(bqd);

   gradrec_type gradrec     = gradrec_type(degree);

   scalar_type error_oper = 0.0;
   scalar_type error_data = 0.0;

   timecounter tc;

   double tgrad = 0.0;
   double tgrad_opt = 0.0;

   for (auto& cl : msh)
   {
      tc.tic();
      gradrec_opt.compute(msh, cl);
      tc.toc();
      tgrad_opt += tc.to_double();

      tc.tic();
      gradrec.compute(msh, cl);
      tc.toc();
      tgrad += tc.to_double();

      error_oper += std::pow( (gradrec_opt.oper() - gradrec.oper).norm()   ,2.0);
      error_data += std::pow( (gradrec_opt.data() - gradrec.data).norm()   ,2.0);

   }

   error_oper = sqrt(error_oper)/msh.cells_size();
   error_data = sqrt(error_data)/msh.cells_size();

   if(error_data > epsilon)
      std::cout << "ecart gradient vector: data " << error_data << std::endl;

   if(error_oper > epsilon)
      std::cout << "ecart gradient vector: oper " << error_oper << std::endl;

   std::cout << "time grad vector: " << tgrad << " vs grad vector opt: " << tgrad_opt << std::endl;

}


template< typename mesh_type, typename T>
void
test_grad_scalar(const mesh_type& msh, const size_t degree, const T epsilon)
{
   typedef typename mesh_type::scalar_type            scalar_type;
   typedef typename mesh_type::cell                   cell_type;
   typedef typename mesh_type::face                   face_type;

   typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
   typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;
   typedef disk::quadrature<mesh_type, cell_type>      matrix_quadrature_type;

   typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
   typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;
   typedef disk::scaled_monomial_vector_basis<mesh_type, cell_type>    matrix_basis_type;

   typedef disk::basis_quadrature_data_full<mesh_type, disk::scaled_monomial_scalar_basis,
   disk::scaled_monomial_vector_basis,
   disk::quadrature> bqdata_type;

   typedef disk::gradient_reconstruction_full_bq<bqdata_type>     gradrecopt_type;


   typedef disk::gradient_reconstruction_full< mesh_type,
                     cell_basis_type, cell_quadrature_type,
                     face_basis_type, face_quadrature_type,
                     matrix_basis_type, matrix_quadrature_type>          gradrec_type;



   typedef dynamic_matrix<scalar_type>         matrix_dynamic;

   bqdata_type bqd(degree, degree, degree);

   gradrecopt_type gradrec_opt(bqd);

   gradrec_type gradrec     = gradrec_type(degree);

   scalar_type error_oper = 0.0;
   scalar_type error_data = 0.0;

   timecounter tc;

   double tgrad = 0.0;
   double tgrad_opt = 0.0;

   for (auto& cl : msh)
   {
      tc.tic();
      gradrec_opt.compute(msh, cl);
      tc.toc();
      tgrad_opt += tc.to_double();

      tc.tic();
      gradrec.compute(msh, cl);
      tc.toc();
      tgrad += tc.to_double();

      error_oper += std::pow( (gradrec_opt.oper() - gradrec.oper).norm()   ,2.0);
      error_data += std::pow( (gradrec_opt.data() - gradrec.data).norm()   ,2.0);

   }

   error_oper = sqrt(error_oper)/msh.cells_size();
   error_data = sqrt(error_data)/msh.cells_size();

   if(error_data > epsilon)
      std::cout << "ecart gradient scalar: data " << error_data << std::endl;

   if(error_oper > epsilon)
      std::cout << "ecart gradient scalar: oper " << error_oper << std::endl;

   std::cout << "time grad scalar: " << tgrad << " vs grad scalar opt: " << tgrad_opt << std::endl;

}




int main(int argc, char **argv)
{
   using RealType = double;

   char    *mesh_filename  = nullptr;
   int     degree          = 1;
   RealType tole(1E-15);

   int ch;

   while ( (ch = getopt(argc, argv, "k:h")) != -1 )
   {
      switch(ch)
      {
         case 'k':
            degree = atoi(optarg);
            if (degree < 0)
            {
               std::cout << "Degree must be positive. Falling back to 1." << std::endl;
               degree = 1;
            }
            break;

         case 'h':
         case '?':
         default:
            std::cout << "wrong arguments" << std::endl;
            exit(1);
      }
   }

   argc -= optind;
   argv += optind;

   if (argc == 0)
   {
      std::cout << "Mesh format: 1D uniform (Not avaible)" << std::endl;
      //auto msh = disk::load_uniform_1d_mesh<RealType>(0, 1, elems_1d);
      //run_NL_elasticity_solver(msh, rp);
      return 0;
   }

   mesh_filename = argv[0];

   /* FVCA5 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
   {
      std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
      auto msh = disk::load_fvca5_2d_mesh<RealType>(mesh_filename);
      test_grad_vector(msh, degree, tole);
      test_grad_scalar(msh, degree, tole);
      return 0;
   }

   /* Netgen 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.mesh2d$") ))
   {
      std::cout << "Guessed mesh format: Netgen 2D" << std::endl;
      auto msh = disk::load_netgen_2d_mesh<RealType>(mesh_filename);
      test_grad_vector(msh, degree, tole);
      test_grad_scalar(msh, degree, tole);
      return 0;
   }

   /* DiSk++ cartesian 2D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
   {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
      auto msh = disk::load_cartesian_2d_mesh<RealType>(mesh_filename);
      test_grad_vector(msh, degree, tole);
      test_grad_scalar(msh, degree, tole);
      return 0;
   }

   /* Netgen 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
   {
      std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
      auto msh = disk::load_netgen_3d_mesh<RealType>(mesh_filename);
      test_grad_vector(msh, degree, tole);
      test_grad_scalar(msh, degree, tole);
      return 0;
   }

   /* DiSk++ cartesian 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.hex$") ))
   {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
      auto msh = disk::load_cartesian_3d_mesh<RealType>(mesh_filename);
      test_grad_vector(msh, degree, tole);
      test_grad_scalar(msh, degree, tole);
      return 0;
   }

   /* DiSk++ cartesian 3D */
   if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
   {
      std::cout << "Guessed mesh format: DiSk++ Cartesian 3D" << std::endl;
      auto msh = disk::load_fvca6_3d_mesh<RealType>(mesh_filename);
      test_grad_vector(msh, degree, tole);
      test_grad_scalar(msh, degree, tole);
      return 0;
   }

   std::cout << "Unkwnon mesh format" << std::endl;
   return 0;

}
