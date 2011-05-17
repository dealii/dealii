//----------------------------  parameter_handler_10.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  parameter_handler_10.cc  ---------------------------


// testing reading a parameter that doesn't conform to its specs. this
// incidentally also uncovered a bug in Patterns::Integer::match

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <fstream>

void check (const char *p)
{
  ParameterHandler prm;
  prm.declare_entry ("test_1", "3", Patterns::Integer());

  std::ifstream in(p);
  const bool result = prm.read_input (in);

  deallog << result << std::endl;
  deallog << "test_1=" << prm.get ("test_1") << std::endl;
}


int main () 
{
  std::ofstream logfile("parameter_handler_10/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check ("parameter_handler_10/prm");
  
  return 0;
}
