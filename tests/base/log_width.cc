// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2014 by the deal.II authors
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


// test that we can set the width of LogStream objects

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>
#include <limits>

int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog << numbers::PI*2 << 42 << std::endl;

  // test with a different width
  deallog.width(18);
  deallog << numbers::PI*2 << 42 << std::endl;

  // ensure that the width of the underlying file stream object remained
  // unchanged
  deallog.get_file_stream() << numbers::PI*2 << 42 << std::endl;

  return 0;
}

