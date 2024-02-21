// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// like _02, but use the parse_input_from_string function

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


int
main()
{
  initlog();

  ParameterHandler prm;
  prm.enter_subsection("Testing");
  prm.declare_entry("Function",
                    "a",
                    Patterns::List(Patterns::Selection("a|b|c|d|e|f|g|h")));
  prm.leave_subsection();

  std::string   input;
  std::ifstream in(SOURCE_DIR "/prm/parameter_handler_2_read_from_string.prm");
  while (in)
    {
      std::string s;
      std::getline(in, s);
      input += s;
      input += '\n';
    }
  prm.parse_input_from_string(input.c_str());

  std::string list;
  prm.enter_subsection("Testing");
  list = prm.get("Function");
  prm.leave_subsection();

  deallog << list << std::endl;

  return 0;
}
