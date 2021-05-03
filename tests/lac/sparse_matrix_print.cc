// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Check that print_as_numpy_arrays works

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"


template <typename number>
void
print_sparse_matrix()
{
  unsigned int const n = 20;
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
