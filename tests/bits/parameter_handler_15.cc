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



// ParameterHandler could not deal with missing endline at end of file
// or can it?
// http://code.google.com/p/dealii/issues/detail?id=126

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <fstream>

void check (const char * content, double &v1, double &v2)
{
  ParameterHandler foo;
  foo.enter_subsection("bar");
  foo.declare_entry("val", "1.0", dealii::Patterns::Double(), "");
  foo.leave_subsection();
  foo.declare_entry("val2", "2.0", dealii::Patterns::Double(), "");

  std::stringstream ss(content);

  foo.read_input(ss);



  foo.enter_subsection("bar");
  deallog << foo.get ("val") << std::endl;
  v1 = foo.get_double("val");
  foo.leave_subsection();
  deallog << foo.get ("val2") << std::endl;
  v2 = foo.get_double("val2");
}

void test(std::string content)
{
  double v1,v2;
  check((content+"\n").c_str(),v1,v2);
  double v3,v4;
  check(content.c_str(),v3,v4);

  Assert(v1==v3, ExcInternalError());
  Assert(v2==v4, ExcInternalError());
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ("subsection bar\nend");
  test ("");
  test ("set val2=-3");
  test ("subsection bar\n set val=2\nend");

  return 0;
}
