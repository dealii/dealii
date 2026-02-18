// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

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
