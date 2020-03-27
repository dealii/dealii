// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Check that, when the correct flag is set, cyclic dependencies in a
// substitution map are resolved in dictionary order, and do not cause an
// infinite loop in the resolve_explicit_dependencies() function.

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

  const bool                        force_cyclic_dependency_resolution = true;
  const SD::types::substitution_map sub_vals_resolved =
    SD::resolve_explicit_dependencies(sub_vals_unresolved,
                                      force_cyclic_dependency_resolution);
  std::cout << "Resolved map:" << std::endl;
  SD::Utilities::print_substitution_map(std::cout, sub_vals_resolved);

  deallog << "OK" << std::endl;
}
