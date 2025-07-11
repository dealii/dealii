// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check ParameterHandler::parse_input_from_json with an undefined parameter

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


int
main()
{
  initlog();

  // default values
  double double1 = 1.234;

  ParameterHandler prm;
  prm.enter_subsection("ss1");
  prm.add_parameter("defined double 1", double1, "doc 3");
  prm.leave_subsection();

  // read from json
  std::ifstream in(SOURCE_DIR "/prm/parameter_handler_read_json_05.json");
  try
    {
      prm.parse_input_from_json(in);
    }
  catch (const ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
      deallog
        << "ParameterHandler threw an exception with the following message:"
        << std::endl;
      e.print_info(deallog.get_file_stream());
    }

  return 0;
}
