// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
