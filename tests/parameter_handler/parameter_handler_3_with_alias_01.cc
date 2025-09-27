// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// like _03, but use an alias to reference two parameters using the same name.
// in the _with_alias_01 testcase, the alias doesn't do anything (we just
// declare it, but the input file still references the original name) whereas
// the _with_alias_02 does it the other way around

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


int
main()
{
  try
    {
      initlog();

      ParameterHandler prm;
      prm.enter_subsection("Testing");
      prm.declare_entry("string list",
                        "a",
                        Patterns::List(Patterns::Selection("a|b|c|d|e|f|g|h")),
                        "docs 1");
      prm.declare_entry("int", "1", Patterns::Integer());
      prm.declare_entry("double", "3.1415926", Patterns::Double(), "docs 3");
      prm.declare_alias("int", "int_alias");
      prm.leave_subsection();

      // read and then write parameters
      prm.parse_input(SOURCE_DIR "/prm/parameter_handler_3.prm");
      prm.print_parameters(deallog.get_file_stream(), ParameterHandler::PRM);
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };

  return 0;
}
