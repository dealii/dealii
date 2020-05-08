// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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


// test LAPACKFullMatrix::solve() for Cholesky factorization

/* MWE for size=3 in Octave:
R = [10,2,3; 2, 20, 6; 3, 6, 90]
x = [2; -7; 1]

> R\x
ans =

   0.268693
  -0.385220
   0.027836
*/

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>

#include "../tests.h"

#include "create_matrix.h"


template <typename NumberType>
void
test()
{
  const unsigned int           size = 3;
  LAPACKFullMatrix<NumberType> M(size);
  M.set_property(LAPACKSupport::symmetric);

  M(0, 0) = 10;
  M(0, 1) = 2;
  M(1, 0) = 2;
  M(0, 2) = 3;
  M(2, 0) = 3;
  M(1, 1) = 20;
  M(1, 2) = 6;
  M(2, 1) = 6;
  M(2, 2) = 90;

  Vector<NumberType> x(size), y(size);
  x[0] = 2;
  x[1] = -7;
  x[2] = 1;

  y = x;
  M.compute_cholesky_factorization();
  M.solve(y);
  y.print(deallog.get_file_stream(), 6, false);
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  test<double>();
}
