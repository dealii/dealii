// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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

// test a map with pair with default and custom pattern
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/point.h>

#include <memory>

#include "../tests.h"

using namespace Patterns::Tools;

int
main()
{
  initlog();

  std::map<unsigned int, std::pair<int, int>> a;
  a[3] = std::make_pair(1, 2);
  a[2] = std::make_pair(3, 4);

  auto b = a;
  auto c = a;

  ParameterHandler prm;
  prm.add_parameter("A map using the default pattern", a);


  prm.add_parameter("A map with explicitly stated default separators",
                    b,
                    "",
                    Patterns::Map(Patterns::Integer(0),
                                  Patterns::Tuple(":",
                                                  Patterns::Integer(),
                                                  Patterns::Integer()),
                                  0,
                                  Patterns::Map::max_int_value,
                                  ",",
                                  "="));

  Patterns::Map custom_pattern(Patterns::Integer(0),
                               Patterns::Tuple(",",
                                               Patterns::Integer(),
                                               Patterns::Integer()),
                               0,
                               Patterns::Map::max_int_value,
                               ";",
                               ":");

  prm.add_parameter("A map with custom separators", c, "", custom_pattern);

  prm.log_parameters(deallog);

  return 0;
}
