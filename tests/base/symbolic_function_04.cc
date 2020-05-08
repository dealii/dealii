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

// Test the SymbolicFunction class: string constructor, multiple components,
// updating constants, and providing additional arguments.

#include <deal.II/base/symbolic_function.h>

#include "../tests.h"

using namespace Differentiation::SD;

int
main()
{
  initlog(0);

  Functions::SymbolicFunction<2> fun("x+beta; alpha*y; t*alpha^2");

  auto sub = make_substitution_map("alpha", Expression("1.0", true));
  fun.update_user_substitution_map(sub);

  auto args = make_substitution_map(make_symbol("beta"), 0.0);
  fun.set_additional_function_arguments(args);

  deallog << fun << std::endl;
  Point<2> p(1, 2);
  deallog << "p: " << p << std::endl
          << "f_0(p): " << fun.value(p, 0) << std::endl
          << "f_1(p): " << fun.value(p, 1) << std::endl
          << "f_2(p): " << fun.value(p, 2) << std::endl;

  sub = make_substitution_map(std::make_pair(make_symbol("alpha"), 2.0));
  fun.update_user_substitution_map(sub);

  deallog << fun << std::endl
          << "p: " << p << std::endl
          << "f_0(p): " << fun.value(p, 0) << std::endl
          << "f_1(p): " << fun.value(p, 1) << std::endl
          << "f_2(p): " << fun.value(p, 2) << std::endl;
}
