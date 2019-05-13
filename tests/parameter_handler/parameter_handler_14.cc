// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2018 by the deal.II authors
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
