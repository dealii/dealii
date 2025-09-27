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



// test ParameterHandler::set(., bool) which was broken, see bug #49

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


int
main()
{
  try
    {
      initlog();

      // same as parameter_handler_3
      ParameterHandler prm;
      prm.enter_subsection("Testing");
      prm.declare_entry("bool", "true", Patterns::Bool(), "docs 1");
      prm.leave_subsection();

      prm.parse_input(SOURCE_DIR "/prm/parameter_handler_6_bool.prm");

      // now set the parameter to a different
      // value
      prm.enter_subsection("Testing");
      prm.set("bool", true);
      prm.leave_subsection();

      // then write
      prm.print_parameters(deallog.get_file_stream(), ParameterHandler::PRM);

      // and do it again with the opposite
      // value
      prm.enter_subsection("Testing");
      prm.set("bool", false);
      prm.leave_subsection();

      // then write
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
