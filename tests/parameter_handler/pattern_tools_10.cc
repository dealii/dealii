// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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

// Check Convert<tuple>::to_value.

#include <deal.II/base/parameter_handler.h>

#include <memory>

#include "../tests.h"

int
main()
{
  initlog();

  auto        a = std::make_tuple(Point<3>(), double(3.5), std::string("ciao"));
  const auto &pattern = Patterns::Tools::Convert<decltype(a)>::to_pattern();

  deallog << pattern->description() << std::endl;

  deallog << Patterns::Tools::Convert<decltype(a)>::to_string(a) << std::endl;

  a = Patterns::Tools::Convert<decltype(a)>::to_value(
    "1.0, 2.0, 3.0 : 3.14 : mondo");

  deallog << Patterns::Tools::Convert<decltype(a)>::to_string(a) << std::endl;
}
