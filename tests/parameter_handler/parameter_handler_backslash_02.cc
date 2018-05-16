// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include "../tests.h"
#include <deal.II/base/parameter_handler.h>

/*
 * Test that ParameterHandler can read parameters of the form
 *
 *     set Function = a,\
 *                    \
 *                    b,\
 *     \
 *                    c
 *
 * correctly. This tests how ParameterHandler handles lines that are blank
 * aside from whitespace and a continuation ('\') character.
 */

int
main ()
{
  initlog();

  for (unsigned int i = 0; i < 2; ++i)
    {
      ParameterHandler prm;
      prm.enter_subsection ("Testing");
      prm.declare_entry ("Function",
                         "a",
                         Patterns::List(Patterns::Selection("a|b|c|d|e|f|g|h")));
      prm.leave_subsection ();

      // test both relevant parse_input functions
      if (i == 0)
        {
          prm.parse_input(SOURCE_DIR "/prm/parameter_handler_backslash_02.prm");
        }
      else
        {
          std::ifstream input_stream
          (SOURCE_DIR "/prm/parameter_handler_backslash_02.prm");
          prm.parse_input(input_stream);
        }

      std::string list;
      prm.enter_subsection ("Testing");
      list = prm.get ("Function");
      prm.leave_subsection ();

      deallog << list << std::endl;
    }

  return 0;
}
