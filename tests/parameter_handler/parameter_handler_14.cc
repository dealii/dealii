// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// ParameterHandler could not deal with parameters named "value" as well as a
// few other names. see the thread on the mailing starting with a post by
// Denis Davydov on March 30, 2013

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

void
check()
{
  std::string input = "subsection bar\n"
                      "  set value = 1.0\n"
                      "end\n";

  ParameterHandler foo;
  foo.enter_subsection("bar");
  foo.declare_entry("value", "1.0", dealii::Patterns::Double(), "");
  foo.leave_subsection();

  try
    {
      foo.parse_input_from_string(input.c_str());
    }
  catch (...)
    {
      deallog << "Exception caught, but none should happen here!!!"
              << std::endl;
    }

  foo.enter_subsection("bar");
  deallog << foo.get("value") << std::endl;
  foo.leave_subsection();
}


int
main()
{
  initlog();

  check();

  return 0;
}
