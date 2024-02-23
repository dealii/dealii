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
