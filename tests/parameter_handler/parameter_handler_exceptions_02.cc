// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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



// ensure that we end up in a defined state if an action throws an
// exception

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


std::string input = "set test_1 = 1\n"
                    "subsection subsec\n"
                    "  set test_2 = -1\n"
                    "end\n";



void
check(const char *p)
{
  ParameterHandler prm;
  prm.declare_entry("test_1", "0", Patterns::Integer(-1, 1));
  prm.enter_subsection("subsec");
  prm.declare_entry("test_2", "0", Patterns::Integer(-1, 1));
  prm.add_action("test_2", [](const std::string &s) {
    // throw an exception from the action for
    // everything but the default value
    if (s != "0")
      throw 1;
  });
  prm.leave_subsection();

  std::istringstream in(input);
  try
    {
      deallog << "Trying to read parameters..." << std::endl;
      prm.parse_input(in);
      deallog << "Done reading parameters..." << std::endl;
    }
  catch (...)
    {
      deallog << "Caught an exception -- ignoring..." << std::endl;
    }


  // make sure the prm object was reset to a state where we are in the
  // subsection we were in before attempting the `parse_input` call
  // (namely, in the top-level section of the prm tree)
  deallog << "test_1="
          << prm.get("test_1") // should =1, because we read that value
          << std::endl;
  prm.enter_subsection("subsec");
  deallog << "test_2="
          << prm.get("test_2") // should =default, because the action failed
          << std::endl;
  prm.leave_subsection();
}


int
main()
{
  initlog();

  check(SOURCE_DIR "/prm/parameter_handler_1.prm");

  return 0;
}
