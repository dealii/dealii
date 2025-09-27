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



// check FullMatrix::swap_col and FullMatrix::swap_row for nonsymmetric
// matrices; there used to be a bug in the implementation


#include "../tests.h"

#include "full_matrix_common.h"



template <typename number>
void
check()
{
  FullMatrix<number> m(5, 7);
  fill_matrix(m);
  print_matrix(m);

  m.swap_col(2, 4);
  print_matrix(m);

  m.swap_row(2, 4);
  print_matrix(m);
}
