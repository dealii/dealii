// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check the Patterns::List pattern with a separator different from the
// default ','

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

void
check(const char *p)
{
  ParameterHandler prm;
  prm.declare_entry("test_1",
                    "-1 xyz 0",
                    Patterns::List(Patterns::Integer(-1, 1), 2, 3, "xyz"));

  std::ifstream in(p);
  prm.parse_input(in);

  deallog << "test_1=" << prm.get("test_1") << std::endl;
}


int
main()
{
  initlog();

  check(SOURCE_DIR "/prm/parameter_handler_1a.prm");

  return 0;
}
