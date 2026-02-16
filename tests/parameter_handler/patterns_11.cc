// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/base/parameter_handler.h>

#include <memory>

#include "../tests.h"

int
main()
{
  initlog();

  // create a pattern and match a string
  const auto &pattern =
    Patterns::Tuple(Patterns::Double(), Patterns::Anything());
  const std::string desc = pattern.description();

  deallog << desc << std::endl;

  std::string test = "3.14: Ciao";

  if (pattern.match(test))
    deallog << "OK" << std::endl;
  else
    deallog << "Not OK" << std::endl;
}
