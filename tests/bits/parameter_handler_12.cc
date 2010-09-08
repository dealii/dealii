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


// test that using ParameterHandler::set with a parameter that doesn't conform
// to the specs leads to an error

#include "../tests.h"
#include <base/logstream.h>
#include <base/parameter_handler.h>
#include <fstream>

void check ()
{
  ParameterHandler prm;
  prm.declare_entry ("test_1", "3", Patterns::Integer());

  try
    {
      prm.set ("test_1", "3.1415");
    }
  catch (const ParameterHandler::ExcValueDoesNotMatchPattern &)
    {
      deallog << "OK" << std::endl;
    }
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
