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



// check creation and output of a matrix using an array


#include "../tests.h"

#include "full_matrix_common.h"



template <typename number>
void
check()
{
  const number array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3,
                          4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5};

  FullMatrix<number> m(5, 5, array);

  print_matrix(m);
}
