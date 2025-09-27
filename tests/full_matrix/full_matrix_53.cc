// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check FullMatrix::left_invert and FullMatrix::right_invert


#include "../tests.h"

#include "full_matrix_common.h"



template <typename number>
void
calculate(const FullMatrix<number> A, const FullMatrix<number> B)
{
  FullMatrix<number> A_r_inv(A.n(), A.m());
  FullMatrix<number> identity(A.m(), A.m());
  A_r_inv.right_invert(A);
  deallog << "A matrix" << std::endl;
  display_matrix(A);
  deallog << "Right inverse" << std::endl;
  display_matrix(A_r_inv);
  deallog << "Identity = A*A_r_inv" << std::endl;
  A.mmult(identity, A_r_inv);
  display_matrix(identity);

  deallog << std::endl;

  FullMatrix<number> B_l_inv(B.n(), B.m());
  FullMatrix<number> identity2(B.n(), B.n());
  B_l_inv.left_invert(B);
  deallog << "B matrix" << std::endl;
  display_matrix(B);
  deallog << "Left inverse" << std::endl;
  display_matrix(B_l_inv);
  deallog << "Identity = B_l_inv*B" << std::endl;
  B_l_inv.mmult(identity2, B);
  display_matrix(identity2);
}


template <typename number>
void
check()
{
  FullMatrix<number> A(1, 2);
  fill_matrix(A);

  FullMatrix<number> B(2, 1);
  fill_matrix(B);

  calculate(A, B);

  FullMatrix<number> A1(2, 3);
  fill_matrix(A1);

  FullMatrix<number> B1(3, 2);
  fill_matrix(B1);

  calculate(A1, B1);
}
