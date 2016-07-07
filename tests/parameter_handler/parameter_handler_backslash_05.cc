// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2016 by the deal.II authors
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
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <fstream>

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

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  for (unsigned int i = 0; i < 2; ++i)
    {
      ParameterHandler prm;
      prm.enter_subsection ("Testing");
      prm.declare_entry ("Function_1",
                         "a",
                         Patterns::List(Patterns::Selection("a|b|c")));
      prm.declare_entry ("Function_2",
                         "d",
                         Patterns::List(Patterns::Selection("d|e|f")));
      prm.leave_subsection ();


      // test both relevant read_input functions
      if (i == 0)
        {
          prm.read_input(SOURCE_DIR "/prm/parameter_handler_backslash_05.prm");
        }
      else
        {
          std::ifstream input_stream
          (SOURCE_DIR "/prm/parameter_handler_backslash_05.prm");
          prm.read_input(input_stream);
        }

      std::string list_1;
      std::string list_2;
      prm.enter_subsection ("Testing");
      list_1 = prm.get ("Function_1");
      list_2 = prm.get ("Function_2");
      prm.leave_subsection ();

      deallog << list_1 << std::endl;
      deallog << list_2 << std::endl;
    }

  return 0;
}
