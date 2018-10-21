// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2018 by the deal.II authors
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


// Tests LAPACKFullMatrix::set(i,j,a)

#include <deal.II/lac/lapack_full_matrix.h>

#include <iostream>
#include <tuple>

#include "../tests.h"

template <typename NumberType>
void
test(const unsigned int n = 3, const unsigned int k = 6)
{
  LAPACKFullMatrix<NumberType> A(n, k);
  A.set(0, 1, 2.);
  DEAL_II_AssertThrow(A(0, 1) == 2., ExcInternalError());

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  test<double>();
}
