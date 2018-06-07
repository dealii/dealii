// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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



// This tests the construction of a BlockMatrixArray and outputs the
// entered blocks using print_latex.

#include <deal.II/lac/block_matrix_array.h>
#include <deal.II/lac/full_matrix.h>

#include "../tests.h"


int
main()
{
  initlog();

  FullMatrix<double> A1(4, 4);
  FullMatrix<double> A2(4, 4);
  FullMatrix<double> B(4, 3);
  FullMatrix<double> C(3, 3);

  BlockMatrixArray<double> block(2, 2);

  block.enter(A1, 0, 0);
  block.enter(A2, 0, 0, 2, true);
  block.enter(B, 0, 1, -3.);
  block.enter(B, 0, 1, -3., true);
  block.enter(C, 1, 1, 1., true);

  block.print_latex(deallog);

  return 0;
}
