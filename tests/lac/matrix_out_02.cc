// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// like matrix_out.cc, but test for Trilinos matrices
//
// also test some of the other options of the MatrixOut::Options class


#include <deal.II/lac/matrix_out.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include "../tests.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);

  // test for a rectangular sparse
  // matrix
  if (true)
    {
      TrilinosWrappers::SparsityPattern sparsity(4, 8, 7);
      for (unsigned int i = 0; i < 4; ++i)
        for (unsigned int j = 0; j < 8; ++j)
          if (i == j + 1)
            sparsity.add(i, j);
      sparsity.compress();

      TrilinosWrappers::SparseMatrix sparse_matrix(sparsity);
      for (unsigned int i = 0; i < 4; ++i)
        for (unsigned int j = 0; j < 8; ++j)
          if (i == j + 1)
            sparse_matrix.set(i, j, i + 3 * j);

      MatrixOut matrix_out;
      matrix_out.build_patches(sparse_matrix,
                               "sparse_matrix",
                               MatrixOut::Options(true, 1, true, false));
      matrix_out.write_gnuplot(logfile);
    }
}
