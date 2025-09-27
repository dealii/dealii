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



// another check FullMatrix::Tadd
// proper transposition of rectangular matrices is verified


#include "../tests.h"

#include "full_matrix_common.h"



template <typename number>
void
check()
{
  FullMatrix<number> A(8, 9);
  fill_matrix(A);
  deallog << "Given matrix" << std::endl;
  display_matrix(A);
  FullMatrix<number> A_t(A.n(), A.m());
  A_t.Tadd(A, 1);
  deallog << "Transposed matrix" << std::endl;
  display_matrix(A_t);
}
