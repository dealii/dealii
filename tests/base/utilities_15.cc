// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test Utilities::type_to_string()

#include <deal.II/base/utilities.h>

#include "../tests.h"

int
main()
{
  initlog();

  double   a = 1.0;
  int      b = 3;
  Point<2> c;

  deallog << Utilities::type_to_string(a) << std::endl
          << Utilities::type_to_string(b) << std::endl
          << Utilities::type_to_string(c) << std::endl;
}
