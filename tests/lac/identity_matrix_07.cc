// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
