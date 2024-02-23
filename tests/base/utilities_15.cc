// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
