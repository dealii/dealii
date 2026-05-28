// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// test LinearAlgebra::TpetraWrappers::SparseMatrix<double>::reinit with another
// LinearAlgebra::TpetraWrappers::SparseMatrix<double>

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_sparsity_pattern.h>

#include <iostream>

#include "../tests.h"


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  initlog();

  SparsityPattern sparsity(5, 5, 5);
  sparsity.add(1, 2);
  sparsity.add(2, 3);
  sparsity.add(3, 4);
  sparsity.add(4, 3);
  sparsity.compress();
  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default>
    matrix;
  matrix.reinit(sparsity);
  {
    double value = 1;
    for (auto p = matrix.begin(); p != matrix.end();
         ++p, ++value)
      p->value() = value;
  }
  matrix.compress(VectorOperation::insert);
  
  deallog << "Original:" << std::endl;
  matrix.print(deallog.get_file_stream());

  // now reinitialize a second Trilinos matrix
  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default>
    tmatrix;
  tmatrix.reinit(matrix);

  deallog << "Copy structure only:" << std::endl;
  tmatrix.print(deallog.get_file_stream());
}
