// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2017 by the deal.II authors
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
