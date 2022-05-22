// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2022 by the deal.II authors
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

  prm.print_parameters(deallog.get_file_stream(), ParameterHandler::Text);
}


int
main()
{
  initlog();
  check();
  return 0;
}
