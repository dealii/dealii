// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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
      deallog << foo.get_double("val") << ' ';
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
