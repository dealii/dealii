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


// test that we can set the width of LogStream objects

#include <limits>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  deallog << numbers::PI * 2 << 42 << std::endl;

  // test with a different width
  deallog.width(18);
  deallog << numbers::PI * 2 << 42 << std::endl;

  // ensure that the width of the underlying file stream object remained
  // unchanged
  deallog.get_file_stream() << numbers::PI * 2 << 42 << std::endl;

  return 0;
}
