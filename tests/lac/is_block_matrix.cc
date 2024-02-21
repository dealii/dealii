// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix_ez.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <algorithm>

#include "../tests.h"



void
test()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);

  deallog
    << internal::AffineConstraints::IsBlockMatrix<SparseMatrix<double>>::value
    << ' '
    << internal::AffineConstraints::IsBlockMatrix<SparseMatrix<float>>::value
    << ' '
    << internal::AffineConstraints::IsBlockMatrix<SparseMatrixEZ<double>>::value
    << ' '
    << internal::AffineConstraints::IsBlockMatrix<SparseMatrixEZ<float>>::value
    << std::endl;

  deallog << internal::AffineConstraints::IsBlockMatrix<
               BlockSparseMatrix<double>>::value
          << ' '
          << internal::AffineConstraints::IsBlockMatrix<
               BlockSparseMatrix<float>>::value
          << ' '
          << internal::AffineConstraints::IsBlockMatrix<
               BlockSparseMatrixEZ<double>>::value
          << ' '
          << internal::AffineConstraints::IsBlockMatrix<
               BlockSparseMatrixEZ<float>>::value
          << std::endl;
}



int
main()
{
  try
    {
      test();
    }
  catch (const std::exception &e)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << e.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      // abort
      return 2;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      // abort
      return 3;
    };


  return 0;
}
