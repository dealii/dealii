// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

// Copying from a matrix into std::vector and vice versa using iterators

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"

int
main()
{
  initlog();

  dealii::FullMatrix<double> A(3, 4);
  std::iota(A.begin(), A.end(), 1.0);

  deallog << "A =" << std::endl;
  A.print(deallog);

  std::vector<double> v(A.n_elements());
  std::copy(A.begin(), A.end(), v.begin());
  deallog << "v = " << std::endl;
  deallog << v << std::endl;

  dealii::FullMatrix<double> B(A.m(), A.n());
  std::copy(v.cbegin(), v.cend(), B.begin());
  deallog << "B =" << std::endl;
  A.print(deallog);

  return 0;
}
