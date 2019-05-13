// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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


// test that we can set the precision of LogStream objects

#include <limits>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  deallog << numbers::PI << std::endl;

  // test with a different precision
  deallog.precision(12);
  deallog << numbers::PI << std::endl;

  // ensure that the precision of the underlying file stream object remained
  // unchanged
  deallog.get_file_stream() << numbers::PI << std::endl;

  return 0;
}
