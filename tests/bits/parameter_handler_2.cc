//----------------------------  parameter_handler_2.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2003, 2004 by the deal.II authors and Brent Bayley
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  parameter_handler_2.cc  ---------------------------


// check the Patterns::List pattern. this particular test failed at
// one point in time with an assertion due to a pretty stupid bug.

#include "../tests.h"
#include <base/logstream.h>
#include <base/parameter_handler.h>
#include <fstream>


int main () 
{
  std::ofstream logfile("parameter_handler_2.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  ParameterHandler prm;
  prm.enter_subsection ("Testing");
  prm.declare_entry ("Function",
                     "a",
                     Patterns::List(Patterns::Selection("a|b|c|d|e|f|g|h")));
  prm.leave_subsection ();

  prm.read_input("parameter_handler_2.prm1");

  std::string list;
  prm.enter_subsection ("Testing");
  list = prm.get ("Function");
  prm.leave_subsection ();

  deallog << list << std::endl;
  
  return 0;
}
