// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that cyclic dependencies in a substitution map are are detected when
// using the resolve_explicit_dependencies() function.

#include <deal.II/differentiation/sd.h>

#include <iostream>
#include <string>

#include "../tests.h"

namespace SD = Differentiation::SD;

int
main()
{
  initlog();

  // Introduce a cyclic dependency into a map.
  // Here we set a==b and b==a.
  const SD::types::substitution_map sub_vals_unresolved =
    SD::make_substitution_map(
      std::make_pair(SD::Expression("a"), SD::Expression("b")),
      std::make_pair(SD::Expression("b"), SD::Expression("a")));
  std::cout << "Original map:" << std::endl;
  SD::Utilities::print_substitution_map(std::cout, sub_vals_unresolved);

  deal_II_exceptions::disable_abort_on_exception();
  try
    {
      const SD::types::substitution_map sub_vals_resolved =
        SD::resolve_explicit_dependencies(sub_vals_unresolved);
      std::cout << "Resolved map:" << std::endl;
      SD::Utilities::print_substitution_map(std::cout, sub_vals_resolved);
      deallog
        << "Cyclic dependency in substitution map was unexpectedly resolved."
        << std::endl;
    }
  catch (...)
    {
      deallog << "Detected cyclic dependency in substitution map." << std::endl;
    }

  deallog << "OK" << std::endl;
}
