// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// testing reading a parameter that doesn't conform to its specs. this
// incidentally also uncovered a bug in Patterns::Integer::match

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

void
check(const char *p)
{
  ParameterHandler prm;
  prm.declare_entry("test_1", "3", Patterns::Integer());

  std::ifstream in(p);
  try
    {
      prm.parse_input(in);

      // The first line in the parameter file should not match the given
      // pattern, so we should not get here
      deallog << "test_1=" << prm.get("test_1") << std::endl;
    }
  catch (ParameterHandler::ExcInvalidEntryForPattern &exc)
    {
      deallog << exc.get_exc_name() << std::endl;
      exc.print_info(deallog.get_file_stream());
    }
}


int
main()
{
  initlog();

  check(SOURCE_DIR "/prm/parameter_handler_10.prm");

  return 0;
}
