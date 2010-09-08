//----------------------------  parameter_handler_12.cc  ---------------------------
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
//----------------------------  parameter_handler_12.cc  ---------------------------


// testing reading a parameter that doesn't conform to its specs. this
// incidentally also uncovered a bug in Patterns::Integer::match

#include "../tests.h"
#include <base/logstream.h>
#include <base/parameter_handler.h>
#include <fstream>

void check ()
{
  ParameterHandler prm;
  prm.declare_entry ("test_1", "3", Patterns::Integer());

  prm.set ("test_1", "3.1415");
  deallog << "test_1=" << prm.get ("test_1") << std::endl;
}


int main () 
{
  std::ofstream logfile("parameter_handler_12/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check ();
  
  return 0;
}
