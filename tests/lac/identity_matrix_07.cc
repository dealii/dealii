// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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


// Check that clearing a FullMatrix<Number> using a number type convertible to
// Number is not ambiguous. Previously, IdentityMatrix(size_type) was detected
// as viable overload.

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/identity_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

int
main()
{
  initlog();
  FullMatrix<std::complex<double>> matrix;
  matrix = 0.;

  deallog << "OK" << std::endl;
  return 0;
}
