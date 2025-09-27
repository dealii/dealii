// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check the Patterns::Map pattern with a separator other than the default ','

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

void
check(const char *p)
{
  ParameterHandler prm;
  prm.declare_entry("test_13",
                    "-1:a xyz 0:b xyz 1:c",
                    Patterns::Map(Patterns::Integer(-1, 1),
                                  Patterns::Selection("a|b|c"),
                                  2,
                                  3,
                                  "xyz"));

  std::ifstream in(p);
  prm.parse_input(in);

  deallog << "test_13=" << prm.get("test_13") << std::endl;
}


int
main()
{
  initlog();

  check(SOURCE_DIR "/prm/parameter_handler_13a.prm");

  return 0;
}
