// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2018 by the deal.II authors
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


// check method mTmult of FullMatrix, symmetric case

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"

const double entries_A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
const double compare[9]   = {14, 32, 50, 32, 77, 122, 50, 122, 194};

int
main()
{
  initlog();
  deallog << std::fixed;
  deallog << std::setprecision(3);

  FullMatrix<double> A(3, 3, entries_A);
  FullMatrix<double> C(3, 3);
  FullMatrix<double> D(3, 3, compare);

  // compute C= A*A^T
  A.mTmult(C, A);

  C.add(-1., D);
  Assert(C.frobenius_norm() < 1e-12, ExcInternalError());

  deallog << "OK" << std::endl;
}
