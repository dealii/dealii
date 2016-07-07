// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// ParameterHandler seemed to ignore everything on the same line behind "end",
// we should generate an error instead.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <fstream>

void check (const char *content)
{
  deallog << "* check" << std::endl;
  ParameterHandler foo;
  foo.enter_subsection("bar");
  foo.declare_entry("val", "1.0", dealii::Patterns::Double(), "");
  foo.leave_subsection();
  foo.declare_entry("val2", "2.0", dealii::Patterns::Double(), "");

  std::stringstream ss(content);

  if (!foo.read_input(ss))
    {
      deallog << "read_input() failed" << std::endl;
      return;
    }

  deallog << "input: ";
  foo.enter_subsection("bar");
  deallog << foo.get_double ("val") << " ";
  foo.leave_subsection();
  deallog << foo.get_double ("val2") << std::endl;
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  check ("subsection bar\nend  #comment is okay");
  check ("subsection bar\nend  ");
  check ("subsection bar\nendhello what is this?");
  check ("subsection bar\nendset val2=-3");
  check ("subsection bar\nendset val2=-3\nset val2=-2");

  return 0;
}
