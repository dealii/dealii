// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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

// Check Convert<tuple>::to_string.

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/std_cxx14/memory.h>

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
}
