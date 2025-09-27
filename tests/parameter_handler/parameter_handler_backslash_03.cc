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

// header for chdir is platform dependent; when the day comes that we support
// Windows in the test suite then conditionally include direction.h
#include <unistd.h>

/*
 * Test that ParameterHandler does not join lines for things like
 *
 *     set Function = a,\ # first term
 *                    b
 *
 * since there are non-whitespace characters after the '\'.
 */

int
main()
{
  initlog();

  for (unsigned int i = 0; i < 2; ++i)
    {
      ParameterHandler prm;
      prm.enter_subsection("Testing");
      prm.declare_entry("Function",
                        "a",
                        Patterns::List(Patterns::Selection("a|b|c|d|e|f|g|h")));
      prm.leave_subsection();

      // We need a local path for the file to get consistent output messages.
      const int chdir_return_code = chdir(SOURCE_DIR);
      AssertThrow(chdir_return_code == 0, ExcInternalError());
      // test both relevant parse_input functions. They should fail with a
      // specific exception.
      try
        {
          if (i == 0)
            {
              prm.parse_input("prm/parameter_handler_backslash_03.prm");
            }
          else
            {
              std::ifstream input_stream(
                "prm/parameter_handler_backslash_03.prm");
              prm.parse_input(input_stream);
            }

          // parse_input should fail and we should not get here
          std::string list;
          prm.enter_subsection("Testing");
          list = prm.get("Function");
          prm.leave_subsection();

          deallog << list << std::endl;
        }
      catch (ParameterHandler::ExcInvalidEntryForPattern &exc)
        {
          deallog << exc.get_exc_name() << std::endl;
          exc.print_info(deallog.get_file_stream());
        }
    }

  return 0;
}
