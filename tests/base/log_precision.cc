// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


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
