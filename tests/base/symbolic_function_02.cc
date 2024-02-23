// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test the SymbolicFunction class: string constructor, and multiple components.

#include <deal.II/base/symbolic_function.h>

#include <map>

#include "../tests.h"

using namespace Differentiation::SD;

int
main()
{
  initlog();

  Functions::SymbolicFunction<2> fun("x; y; t");
  deallog << fun << std::endl;
  Point<2> p(1, 2);
  deallog << "p: " << p << std::endl
          << "f_0(p): " << fun.value(p, 0) << std::endl
          << "f_1(p): " << fun.value(p, 1) << std::endl
          << "f_2(p): " << fun.value(p, 2) << std::endl;
}
