// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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



// ParameterHandler could not deal with parameters named "value" as well as a
// few other names. see the thread on the mailing starting with a post by
// Denis Davydov on March 30, 2013

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <fstream>

void check ()
{
  ParameterHandler foo;
  foo.enter_subsection("bar");
  foo.declare_entry("value", "1.0", dealii::Patterns::Double(), "");
  foo.leave_subsection();

  try
    {
      foo.read_input("tmp.prm");
    }
  catch (...)
    {
      deallog << "Exception caught, but none should happen here!!!"
              << std::endl;
    }

  foo.enter_subsection("bar");
  deallog << foo.get ("value") << std::endl;
  foo.leave_subsection();

  // delete tmp file again
  std::remove("tmp.prm");
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check ();

  return 0;
}
