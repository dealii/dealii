// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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
