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



// check that setting last_line parameter does indeed terminate parsing.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <fstream>

void check (const char *p, std::string last_line)
{
  ParameterHandler prm;
  prm.declare_entry ("var_1", "-1",
                     Patterns::Integer(-1,1));

  prm.declare_entry ("var_2", "-1",
                     Patterns::Integer(-1,1));

  std::ifstream in(p);
  bool status = prm.read_input (in, "input file", last_line);
  Assert (status == true, ExcInternalError());

  deallog << "var_1=" << prm.get ("var_1") << std::endl
          << "var_2=" << prm.get ("var_2") << std::endl;
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  // go into the source dir to read files there. this
  // is necessary so that we can include files there
  chdir (SOURCE_DIR);
  check ("parameter_handler_22_last_line.prm","");
  check ("parameter_handler_22_last_line.prm","# end-of-parameter-section");

  return 0;
}
