// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2006 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// similar to log_nan, but test for infinities

#include <limits>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  deallog << std::numeric_limits<double>::infinity() << std::endl;
  deallog << -std::numeric_limits<double>::infinity() << std::endl;

  return 0;
}
