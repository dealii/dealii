// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// test LAPACKFullMatrix removing row and column from a matrix

#include "../tests.h"
#include <deal.II/lac/lapack_full_matrix.h>

#include <iostream>


template <typename NumberType>
void
test()
{
  const unsigned int n_rows = 6;
  const unsigned int n_cols = 9;
  LAPACKFullMatrix<NumberType> A(n_rows,n_cols);

  for (unsigned int i = 0; i < n_rows; ++i)
    for (unsigned int j = 0; j < n_cols; ++j)
      A(i,j) = 10*(i+1)+(j+1);

  A.remove_row_and_column(2,7); // <-- 3x and x8 are removed
  A.print_formatted(deallog.get_file_stream(), 0, false, 2);
}


int
main()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  logfile.precision(5);
  deallog.attach(logfile);

  test<double>();
}
