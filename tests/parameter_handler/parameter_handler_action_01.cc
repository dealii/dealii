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



// check the ParameterHandler::add_action() function


#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


void
check(const char *p)
{
  std::string parameter_set_by_action;

  ParameterHandler prm;
  prm.declare_entry("test_1",
                    "-1,0",
                    Patterns::List(Patterns::Integer(-1, 1), 2, 3));
  prm.add_action("test_1", [&](const std::string &s) {
    deallog << "In action:" << s << std::endl;
    parameter_set_by_action = s;
    return true;
  });


  std::ifstream in(p);

  deallog << "Reading parameters" << std::endl;
  prm.parse_input(in);

  deallog << "test_1=" << prm.get("test_1") << std::endl;
  deallog << "Saved parameter: " << parameter_set_by_action << std::endl;
}



int
main()
{
  initlog();

  check(SOURCE_DIR "/prm/parameter_handler_1.prm");

  return 0;
}
