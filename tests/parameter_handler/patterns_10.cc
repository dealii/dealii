// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/parameter_handler.h>

#include <memory>

#include "../tests.h"

int
main()
{
  initlog();

  // create a pattern and match a string
  const auto       &pattern = Patterns::Tuple(";",
                                        Patterns::Integer(),
                                        Patterns::Double(),
                                        Patterns::Anything());
  const std::string desc    = pattern.description();

  deallog << desc << std::endl;

  std::string test = "5; 3.14; Ciao";

  if (pattern.match(test))
    deallog << "OK" << std::endl;
  else
    deallog << "Not OK" << std::endl;
}
