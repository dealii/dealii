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



// check that we can do include statements. the current test verifies what
// happens if such an include statement fails

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

void
check(const char *p)
{
  ParameterHandler prm;
  prm.declare_entry("test_1",
                    "-1,0",
                    Patterns::List(Patterns::Integer(-1, 1), 2, 3));

  std::ifstream in(p);
  try
    {
      prm.parse_input(in);
    }
  catch (ParameterHandler::ExcCannotOpenIncludeStatementFile &exc)
    {
      deallog << exc.get_exc_name() << std::endl;
      exc.print_info(deallog.get_file_stream());
    }

  // Even though the parameter handler failed to finish parsing, it should
  // have still picked up the first statement:
  deallog << "test_1=" << prm.get("test_1") << std::endl;
}


int
main()
{
  initlog();

  // go into the source dir to read files there. this
  // is necessary so that we can include files there
  chdir(SOURCE_DIR);
  check("parameter_handler_1_include_fail.prm");

  return 0;
}
