// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// ensure that we end up in a defined state after a pattern is not matched

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


std::string input = "set test_1 = 1\n"
                    "subsection subsec\n"
                    "  set test_2 = 42\n" // forbidden
                    "end\n";



void
check(const char *p)
{
  ParameterHandler prm;
  prm.declare_entry("test_1", "0", Patterns::Integer(-1, 1));
  prm.enter_subsection("subsec");
  prm.declare_entry("test_2", "0", Patterns::Integer(-1, 1));
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
          << prm.get(
               "test_2") // should =default, because reading the value failed
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
