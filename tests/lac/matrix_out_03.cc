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


// Compare dense and sparse visualizations of a sparse matrix.

#include <deal.II/lac/matrix_out.h>
#include <deal.II/lac/sparse_matrix.h>

#include "../tests.h"

int
main()
{
  initlog();

  auto &logfile = deallog.get_file_stream();

  // Create a rectangular sparse matrix
  SparsityPattern sparsity(14, 18, 7);
  for (unsigned int i = 0; i < 14; ++i)
    for (unsigned int j = 0; j < 18; ++j)
      if ((i + j) % 6 == 0)
        sparsity.add(i, j);
  sparsity.compress();

  SparseMatrix<double> sparse_matrix(sparsity);
  for (unsigned int i = 0; i < 14; ++i)
    for (unsigned int j = 0; j < 18; ++j)
      if ((i + j) % 6 == 0)
        sparse_matrix.set(i, j, i + j);

  // Output it as a dense matrix
  {
    MatrixOut matrix_out;
    matrix_out.build_patches(sparse_matrix,
                             "sparse_matrix",
                             MatrixOut::Options(true, 1, true, false));
    matrix_out.write_vtk(logfile);
  }

  // Output it as a dense matrix
  {
    MatrixOut matrix_out;
    matrix_out.build_patches(sparse_matrix,
                             "sparse_matrix",
                             MatrixOut::Options(true, 1, true, true));
    matrix_out.write_vtk(logfile);
  }
}
