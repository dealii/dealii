// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check skip_undefined settings in parameter handler.

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

void
check()
{
  ParameterHandler prm;
  prm.enter_subsection("Geometry");
  prm.declare_entry("dim", "1", Patterns::Integer(1, 3));
  prm.leave_subsection();

  prm.declare_entry("Dimension", "1", Patterns::Integer(1, 3));

  std::ifstream in(SOURCE_DIR "/parameter_handler_24_in.prm");
  prm.parse_input(in, "input file", "", true);

  prm.print_parameters(deallog.get_file_stream(), ParameterHandler::PRM);
}


int
main()
{
  initlog();
  check();
  return 0;
}
