
// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check DynamicSparsityPattern::compute_mmult_pattern(). Test if multiplication
// of two patterns yield the right DynamicSparsityPattern.

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  // create two different FullMatrix objects and add entries randomly to it.
  // Multiply them via FullMatrix::mmult(..), go through each row and column,
  // and check if DynamicSparsityPattern::compute_mmult_pattern(..) predicts
  // the entries (i,j) where values supposed to be because of the
  // multiplication.

  const unsigned int M = 100;
  const unsigned int N = 50;
  const unsigned int O = 75;

  const unsigned int m_entries = 4;
  const unsigned int n_entries = 5;
  const unsigned int o_entries = 6;

  bool test_failed = false;

  DynamicSparsityPattern dyn_left(M, N), dyn_right(N, O);
  DynamicSparsityPattern dyn_out;

  FullMatrix<double> mat_left(M, N), mat_right(N, O), mat_out(M, O);

  // add randomly entries to the left matrices/pattern starting at the 2nd row
  for (unsigned int m = 1; m < M; ++m)
    for (unsigned int n = 0; n < N; ++n)
      if (Utilities::generate_normal_random_number(0, 0.2) > 0)
        {
          dyn_left.add(m, n);
          mat_left[m][n] = 1;
        }
  // add randomly entries to the right matrices/pattern starting at the 2nd
  // column
  for (unsigned int n = 0; n < N; ++n)
    for (unsigned int o = 1; o < O; ++o)
      if (Utilities::generate_normal_random_number(0, 0.2) > 0)
        {
          dyn_right.add(n, o);
          mat_right[n][o] = 1;
        }
  // first row of left and first column of right in such a way that no entry
  // is created in final matrix
  for (unsigned int n = 0; n < N; ++n)
    if (Utilities::generate_normal_random_number(0, 0.2) > 0)
      {
        dyn_right.add(n, 0);
        mat_right[n][0] = 1;
      }
    else
      {
        dyn_left.add(0, n);
        mat_left[0][n] = 1;
      }


  // do matrix-matrix-multiplication
  mat_left.mmult(mat_out, mat_right);
  dyn_out.compute_mmult_pattern(dyn_left, dyn_right);

  // go through whole matrix and check if pattern are equal
  for (unsigned int i = 0; i < M; ++i)
    for (unsigned int j = 0; j < O; ++j)
      {
        if (mat_out[i][j] != 0 && !dyn_out.exists(i, j))
          test_failed = true;
        else if (mat_out[i][j] == 0 && dyn_out.exists(i, j))
          test_failed = true;
      }

  Assert(!test_failed, ExcInternalError());
  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
  return 0;
}
