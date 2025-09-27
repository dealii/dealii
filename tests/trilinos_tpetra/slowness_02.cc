// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this is part of a whole suite of tests that checks the relative speed of
// using Trilinos for sparse matrices as compared to the speed of our own
// library. the tests therefore may not all actually use Trilinos, but they are
// meant to compare it
//
// the tests build the 5-point stencil matrix for a uniform grid of size N*N

#include <deal.II/base/utilities.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>

#include <iostream>

#include "../tests.h"


void
test()
{
  const unsigned int N = 200;

  // build the sparse matrix
  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default>
    matrix(N * N, N * N, 5U);
  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < N; ++j)
      {
        const unsigned int global = i * N + j;
        matrix.set(global, global, 4);
        if (j > 0)
          {
            matrix.set(global - 1, global, -1);
            matrix.set(global, global - 1, -1);
          }
        if (j < N - 1)
          {
            matrix.set(global + 1, global, -1);
            matrix.set(global, global + 1, -1);
          }
        if (i > 0)
          {
            matrix.set(global - N, global, -1);
            matrix.set(global, global - N, -1);
          }
        if (i < N - 1)
          {
            matrix.set(global + N, global, -1);
            matrix.set(global, global + N, -1);
          }
      }
  matrix.compress(VectorOperation::insert);


  // then do a single matrix-vector
  // multiplication with subsequent formation
  // of the matrix norm
  LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> v1;
  v1.reinit(complete_index_set(N * N), MPI_COMM_WORLD);
  LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> v2;
  v2.reinit(complete_index_set(N * N), MPI_COMM_WORLD);
  for (unsigned int i = 0; i < N * N; ++i)
    v1(i) = i;
  matrix.vmult(v2, v1);

  deallog << v1 * v2 << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  try
    {
      {
        test();
      }
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
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
      return 1;
    };
  return 0;
}
