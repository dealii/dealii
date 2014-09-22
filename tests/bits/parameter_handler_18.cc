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

// this is a variant of parameter_handler_15 but in fact reads data
// from a file

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <fstream>

void test ()
{
  ParameterHandler foo;
  foo.enter_subsection("bar");
  foo.declare_entry("val", "1.0", dealii::Patterns::Double(), "");
  foo.leave_subsection();

  bool okay = foo.read_input(SOURCE_DIR "/parameter_handler_18.prm");
  Assert(okay, ExcMessage("read_input failed"));
  
  foo.enter_subsection("bar");
  deallog << foo.get ("val") << std::endl;
  foo.leave_subsection();
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  return 0;
}
