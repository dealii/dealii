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



// ParameterHandler seemed to ignore everything on the same line behind "end",
// we should generate an error instead.

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

void
check(const char *content)
{
  deallog << "* check" << std::endl;
  ParameterHandler foo;
  foo.enter_subsection("bar");
  foo.declare_entry("val", "1.0", dealii::Patterns::Double(), "");
  foo.leave_subsection();
  foo.declare_entry("val2", "2.0", dealii::Patterns::Double(), "");

  std::stringstream ss(content);

  try
    {
      foo.parse_input(ss);
      deallog << "input: ";
      foo.enter_subsection("bar");
      deallog << foo.get_double("val") << " ";
      foo.leave_subsection();
      deallog << foo.get_double("val2") << std::endl;
    }

  catch (ParameterHandler::ExcCannotParseLine &)
    {
      deallog << "parse_input() failed" << std::endl;
    }
}

int
main()
{
  initlog();

  check("subsection bar\nend  #comment is okay");
  check("subsection bar\nend  ");
  check("subsection bar\nendhello what is this?");
  check("subsection bar\nendset val2=-3");
  check("subsection bar\nendset val2=-3\nset val2=-2");

  return 0;
}
