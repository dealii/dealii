// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Check LAPACKFullMatrix::print_formatted on NaN entry

#include <deal.II/lac/lapack_full_matrix.h>

#include "../tests.h"


int
main()
{
  initlog();

  LAPACKFullMatrix<double> matrix(2, 2);
  matrix(0, 0) = std::numeric_limits<double>::quiet_NaN();
  matrix(0, 1) = std::numeric_limits<double>::infinity();
  matrix(1, 1) = -std::numeric_limits<double>::infinity();

  deallog << "Using print_formatted" << std::endl;
  matrix.print_formatted(deallog.get_file_stream(), 3, true, 0, "0");

  return 0;
}
