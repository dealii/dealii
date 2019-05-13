// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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
