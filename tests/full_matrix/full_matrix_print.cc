// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
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


// verify that FullMatrix::print saves the stream state around the output

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"

const double entries[9] =
  {11.1, 12.2, 13.3, 21.456, 22.12345678901, 23, 31, 32, 33};


int
main()
{
  initlog();
  std::ostream &logfile = deallog.get_file_stream();
  deallog << std::fixed;
  deallog << std::setprecision(3);
  logfile << numbers::PI << std::endl;

  FullMatrix<double> T(3, 3, entries);

  // try writing to a regular stream
  T.print(logfile, 15, 8);

  // print something normal -- should use the old settings, 6 digits of
  // precision
  logfile << numbers::PI << std::endl;

  // now try the same with a LogStream
  T.print(deallog, 15, 8);

  logfile << numbers::PI << std::endl;
}
