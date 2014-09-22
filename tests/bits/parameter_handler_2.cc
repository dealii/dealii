// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// check the Patterns::List pattern. this particular test failed at
// one point in time with an assertion due to a pretty stupid bug.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <fstream>


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  ParameterHandler prm;
  prm.enter_subsection ("Testing");
  prm.declare_entry ("Function",
                     "a",
                     Patterns::List(Patterns::Selection("a|b|c|d|e|f|g|h")));
  prm.leave_subsection ();

  prm.read_input(SOURCE_DIR "/prm/parameter_handler_2.prm");

  std::string list;
  prm.enter_subsection ("Testing");
  list = prm.get ("Function");
  prm.leave_subsection ();

  deallog << list << std::endl;

  return 0;
}
