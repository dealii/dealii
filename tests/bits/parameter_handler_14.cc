//----------------------------  parameter_handler_14.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2010, 2012, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  parameter_handler_14.cc  ---------------------------


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
      foo.read_input("parameter_handler_14/tmp.prm");
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
  std::remove("parameter_handler_14/tmp.prm");
}


int main ()
{
  std::ofstream logfile("parameter_handler_14/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check ();

  return 0;
}
