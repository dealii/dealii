// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2022 by the deal.II authors
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

// Check Patterns::Tools::Convert for std::array types

#include <deal.II/base/parameter_handler.h>

#include <memory>

#include "../tests.h"

int
main()
{
  initlog();

  std::array<Point<2>, 2> points{{Point<2>(0, 1), Point<2>(2, 3)}};

  auto s = Patterns::Tools::to_string(points);

  Patterns::Tools::to_value("4,5 ; 6,7", points);

  deallog << "From: " << s << " to " << Patterns::Tools::to_string(points)
          << std::endl;
}
