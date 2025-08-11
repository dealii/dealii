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



#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/matrix_out.h>

#include "../tests.h"

int
main()
{
  initlog();

  auto &logfile = deallog.get_file_stream();

  // test for a square full matrix
  if (true)
    {
      FullMatrix<double> full_matrix(4, 4);
      for (unsigned int i = 0; i < 4; ++i)
        full_matrix(i, i) = 1;

      MatrixOut matrix_out;
      matrix_out.build_patches(full_matrix,
                               "full_matrix",
                               {false, 1, false, false});
      matrix_out.write_gnuplot(logfile);
    };

  // test for a rectangular sparse
  // matrix
  if (true)
    {
      SparsityPattern sparsity(4, 8, 7);
      for (unsigned int i = 0; i < 4; ++i)
        for (unsigned int j = 0; j < 8; ++j)
          if (i != j)
            sparsity.add(i, j);
      sparsity.compress();

      SparseMatrix<double> sparse_matrix(sparsity);
      for (unsigned int i = 0; i < 4; ++i)
        for (unsigned int j = 0; j < 8; ++j)
          sparse_matrix.set(i, j, static_cast<signed int>(i - j));

      MatrixOut matrix_out;
      matrix_out.build_patches(sparse_matrix,
                               "sparse_matrix",
                               MatrixOut::Options(true, 1, false, false));
      matrix_out.write_eps(logfile);
    };

  // test collation of elements
  if (true)
    {
      FullMatrix<double> full_matrix(20, 20);
      for (unsigned int i = 0; i < 20; ++i)
        for (unsigned int j = 0; j < 20; ++j)
          full_matrix(i, j) =
            (1. * i * i / 20 / 20 - 1. * j * j * j / 20 / 20 / 20);

      MatrixOut matrix_out;
      matrix_out.build_patches(full_matrix,
                               "collated_matrix",
                               MatrixOut::Options(false, 4, false, false));
      matrix_out.write_gmv(logfile);
    };
}
