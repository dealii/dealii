// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check creation and output of a matrix


#include "../tests.h"

#include "full_matrix_common.h"



template <typename number>
void
check()
{
  // create 4x4 matrix and fill with consecutive integers
  FullMatrix<number> mat(4, 4);
  number             val = 1;
  for (unsigned int i = 0; i < mat.m(); ++i)
    for (unsigned int j = 0; j < mat.n(); ++j)
      mat(i, j) = val++;

  // define row and column permutations
  std::vector<unsigned int> row_perm = {2,
                                        0,
                                        3,
                                        1}; // new row 0 <- old row 2, ...
  std::vector<unsigned int> col_perm = {1,
                                        3,
                                        0,
                                        2}; // new col 0 <- old col 1, ...

  // print before permutation
  deallog << "Before permute:" << std::endl;
  print_matrix(mat);

  // apply permutation in-place
  mat.permute(row_perm, col_perm);

  // print after permutation
  deallog << "After permute:" << std::endl;
  print_matrix(mat);
}
