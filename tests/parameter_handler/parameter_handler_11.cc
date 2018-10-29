// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2017 by the deal.II authors
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



// like _10 but for Patterns::Double::match : the first line of the parameter
// file does not match the given pattern.

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

void
check(const char *p)
{
  ParameterHandler prm;
  prm.declare_entry("test_1", "3", Patterns::Double());

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

  check(SOURCE_DIR "/prm/parameter_handler_11.prm");

  return 0;
}
