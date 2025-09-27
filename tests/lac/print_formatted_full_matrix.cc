// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check FullMatrix::print_formatted on NaN entry

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"


int
main()
{
  initlog();

  FullMatrix<double> matrix(2, 2);
  matrix(0, 0) = std::numeric_limits<double>::quiet_NaN();
  matrix(0, 1) = std::numeric_limits<double>::infinity();
  matrix(1, 1) = -std::numeric_limits<double>::infinity();

  deallog << "Using print" << std::endl;
  matrix.print(deallog.get_file_stream());

  deallog << "Using print_formatted" << std::endl;
  matrix.print_formatted(deallog.get_file_stream(), 3, true, 0, "0");

  return 0;
}
