// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that FullMatrix objects can be move constructed and assigned

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"

int
main()
{
  initlog();

  std::size_t        m = 2, n = 3;
  FullMatrix<double> A(m, n);
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t j = 0; j < n; ++j)
      A(i, j) = n * i + j;

  deallog << "Size of A:" << std::endl << A.m() << ' ' << A.n() << std::endl;

  FullMatrix<double> B = std::move(A);

  deallog << "Size of B:" << std::endl << B.m() << ' ' << B.n() << std::endl;
  deallog << "Size of A:" << std::endl << A.m() << ' ' << A.n() << std::endl;

  A = std::move(B);
  deallog << "Size of B:" << std::endl << B.m() << ' ' << B.n() << std::endl;
  deallog << "Size of A:" << std::endl << A.m() << ' ' << A.n() << std::endl;

  return 0;
}
