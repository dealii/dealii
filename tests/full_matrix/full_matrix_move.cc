// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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


// Check that FullMatrix objects can be move constructed and assigned

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"

int
main()
{
  initlog();

  size_t             m = 2, n = 3;
  FullMatrix<double> A(m, n);
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j)
      A(i, j) = n * i + j;

  deallog << "Size of A:" << std::endl << A.m() << " " << A.n() << std::endl;

  FullMatrix<double> B = std::move(A);

  deallog << "Size of B:" << std::endl << B.m() << " " << B.n() << std::endl;
  deallog << "Size of A:" << std::endl << A.m() << " " << A.n() << std::endl;

  A = std::move(B);
  deallog << "Size of B:" << std::endl << B.m() << " " << B.n() << std::endl;
  deallog << "Size of A:" << std::endl << A.m() << " " << A.n() << std::endl;

  return 0;
}
