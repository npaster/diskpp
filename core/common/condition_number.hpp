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

#include <iostream>
#include <sstream>

#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/DenseSymShiftSolve.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymEigsShiftSolver.h>

#include "common/eigen.hpp"

// Compute the condition number of real symetric matrix

template<typename T>
T condition_number(const Eigen::SparseMatrix<T>& mat)
{
   Spectra::SparseSymMatProd<T> op_largest(mat);
   Spectra::SparseSymShiftSolve<T> op_smallest(mat);

   // Construct solver object, requesting the largest eigenvalues
   Spectra::SymEigsSolver< T, Spectra::LARGEST_MAGN, Spectra::SparseSymMatProd<T> >
   eigs_largest(&op_largest, 1, 25);

   // Initialize and compute
   eigs_largest.init();
   eigs_largest.compute();

   if(eigs_largest.info() != Spectra::SUCCESSFUL)
   {
      throw std::invalid_argument("Fail to compute largest eigenvalue");
   }
   const auto eig_largest = eigs_largest.eigenvalues();

   // Construct solver object, requesting the smallest eigenvalues
   Spectra::SymEigsShiftSolver< T, Spectra::LARGEST_MAGN, Spectra::SparseSymShiftSolve<T> >
   eigs_smallest(&op_smallest, 1, 25, T(0));

   eigs_smallest.init();
   eigs_smallest.compute();

   if(eigs_smallest.info() != Spectra::SUCCESSFUL)
   {
      throw std::invalid_argument("Fail to compute smallest eigenvalue");
   }
   const auto eig_smallest = eigs_smallest.eigenvalues();

   std::cout << "largest: " << eig_largest(0) << " smallest: " << eig_smallest(0) << std::endl;

   return std::abs(eig_largest(0)/ eig_smallest(0));
}



template<typename T>
T condition_number(const dynamic_matrix<T>& mat)
{
   Spectra::DenseSymMatProd<T> op_largest(mat);
   Spectra::DenseSymShiftSolve<T> op_smallest(mat);

   // Construct solver object, requesting the largest eigenvalues
   Spectra::SymEigsSolver< T, Spectra::LARGEST_MAGN, Spectra::DenseSymMatProd<T> >
   eigs_largest(&op_largest, 1, 25);

   // Initialize and compute
   eigs_largest.init();
   eigs_largest.compute();

   if(eigs_largest.info() != Spectra::SUCCESSFUL)
   {
      throw std::invalid_argument("Fail to compute largest eigenvalue");
   }
   const auto eig_largest = eigs_largest.eigenvalues();

   // Construct solver object, requesting the smallest eigenvalues
   Spectra::SymEigsShiftSolver< T, Spectra::LARGEST_MAGN, Spectra::DenseSymShiftSolve<T> >
   eigs_smallest(&op_smallest, 1, 25, T(0));

   eigs_smallest.init();
   eigs_smallest.compute();

   if(eigs_smallest.info() != Spectra::SUCCESSFUL)
   {
      throw std::invalid_argument("Fail to compute smallest eigenvalue");
   }
   const auto eig_smallest = eigs_smallest.eigenvalues();

   return std::abs(eig_largest(0)/ eig_smallest(0));
}

template<typename T, size_t N>
T condition_number(const static_matrix<T,N,N>& mat)
{
   Spectra::DenseSymMatProd<T> op_largest(mat);
   Spectra::DenseSymShiftSolve<T> op_smallest(mat);

   // Construct solver object, requesting the largest eigenvalues
   Spectra::SymEigsSolver< T, Spectra::LARGEST_MAGN, Spectra::DenseSymMatProd<T> >
   eigs_largest(&op_largest, 1, 25);

   // Initialize and compute
   eigs_largest.init();
   eigs_largest.compute();

   if(eigs_largest.info() != Spectra::SUCCESSFUL)
   {
      throw std::invalid_argument("Fail to compute largest eigenvalue");
   }
   const auto eig_largest = eigs_largest.eigenvalues();

   // Construct solver object, requesting the smallest eigenvalues
   Spectra::SymEigsShiftSolver< T, Spectra::LARGEST_MAGN, Spectra::DenseSymShiftSolve<T> >
   eigs_smallest(&op_smallest, 1, 25, T(0));

   eigs_smallest.init();
   eigs_smallest.compute();

   if(eigs_smallest.info() != Spectra::SUCCESSFUL)
   {
      throw std::invalid_argument("Fail to compute smallest eigenvalue");
   }
   const auto eig_smallest = eigs_smallest.eigenvalues();

   return std::abs(eig_largest(0)/ eig_smallest(0));
}


template<typename T>
T largest_eigenvalue(const Eigen::SparseMatrix<T>& mat)
{
   Spectra::SparseSymMatProd<T> op_largest(mat);
   Spectra::SparseSymShiftSolve<T> op_smallest(mat);

   // Construct solver object, requesting the largest eigenvalues
   Spectra::SymEigsSolver< T, Spectra::LARGEST_MAGN, Spectra::SparseSymMatProd<T> >
   eigs_largest(&op_largest, 1, 10);

   // Initialize and compute
   eigs_largest.init();
   eigs_largest.compute();

   if(eigs_largest.info() != Spectra::SUCCESSFUL)
   {
      throw std::invalid_argument("Fail to compute largest eigenvalue");
   }
   const auto eig_largest = eigs_largest.eigenvalues();
   
   return eig_largest(0);
}