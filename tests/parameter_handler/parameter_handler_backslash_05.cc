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


#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

/*
 * Test that ParameterHandler will stop a line continuation if a completely
 * blank line follows one with a '\', such as
 *
 *     set Function_1 = a, \
 *
 *                      b, \
 *                      c
 *
 * This should *not* be parsed as 'Function_1 = a, b, c'.
 */

int
main()
{
  initlog();

  for (unsigned int i = 0; i < 2; ++i)
    {
      ParameterHandler prm;
      prm.enter_subsection("Testing");
      prm.declare_entry("Function_1",
                        "a",
                        Patterns::List(Patterns::Selection("a|b|c")));
      prm.declare_entry("Function_2",
                        "d",
                        Patterns::List(Patterns::Selection("d|e|f")));
      prm.leave_subsection();


      // We need a local path for the file to get consistent output messages.
      const int chdir_return_code = chdir(SOURCE_DIR);
      AssertThrow(chdir_return_code == 0, ExcInternalError());
      // test both relevant parse_input functions
      try
        {
          if (i == 0)
            {
              prm.parse_input("prm/parameter_handler_backslash_05.prm");
            }
          else
            {
              std::ifstream input_stream(
                "prm/parameter_handler_backslash_05.prm");
              prm.parse_input(input_stream);
            }

          std::string list_1;
          std::string list_2;
          prm.enter_subsection("Testing");
          list_1 = prm.get("Function_1");
          list_2 = prm.get("Function_2");
          prm.leave_subsection();

          deallog << list_1 << std::endl;
          deallog << list_2 << std::endl;
        }
      catch (ParameterHandler::ExcCannotParseLine &exc)
        {
          deallog << exc.get_exc_name() << std::endl;
          exc.print_info(deallog.get_file_stream());
        }
    }

  return 0;
}
