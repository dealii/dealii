// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check method FullMatrix::scatter_matrix_to

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"

void
test()
{
  // create a matrix with known
  // elements
  FullMatrix<double> A(5, 6);
  for (unsigned int i = 0; i < A.m(); ++i)
    for (unsigned int j = 0; j < A.n(); ++j)
      A(i, j) = i + j;

  // pick every other row and column
  std::vector<types::global_dof_index> rows(A.m());
  for (unsigned int i = 0; i < rows.size(); ++i)
    rows[i] = 2 * i;

  std::vector<types::global_dof_index> cols(A.n());
  for (unsigned int i = 0; i < cols.size(); ++i)
    cols[i] = 2 * i;

  // do the scatter
  FullMatrix<double> X(rows.size() * 2, cols.size() * 2);
  A.scatter_matrix_to(rows, cols, X);

  // verify that the elements are
  // correct
  for (unsigned int i = 0; i < X.m(); ++i)
    for (unsigned int j = 0; j < X.n(); ++j)
      {
        if ((i % 2 == 0) && (j % 2 == 0))
          {
            AssertThrow(X(i, j) == i / 2 + j / 2, ExcInternalError());
          }
        else
          {
            AssertThrow(X(i, j) == 0, ExcInternalError());
          }
      }

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test();
}
