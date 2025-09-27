// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that print_as_numpy_arrays works

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"


template <typename number>
void
print_sparse_matrix()
{
  const unsigned int n = 20;
  SparsityPattern    sp(n, n);
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
      sp.add(i, j);

  sp.compress();

  SparseMatrix<double> A(sp);

  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
      A.set(i, j, static_cast<number>(i + j) / 2.);
  A.print_as_numpy_arrays(deallog.get_file_stream());
}


int
main()
{
  initlog();

  print_sparse_matrix<float>();
  print_sparse_matrix<double>();
}
