// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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
  const int chdir_return_code = chdir(SOURCE_DIR);
  AssertThrow(chdir_return_code == 0, ExcInternalError());
  check("parameter_handler_1_include_fail_in.prm");

  return 0;
}
