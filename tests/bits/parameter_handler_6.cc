//----------------------------  parameter_handler_3.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  parameter_handler_3.cc  ---------------------------


// test ParameterHandler::set(Text)

#include "../tests.h"
#include <base/logstream.h>
#include <base/parameter_handler.h>
#include <fstream>
#include <iomanip>


int main () 
{
  try 
    {
      std::ofstream logfile("parameter_handler_6/output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

                                       // same as parameter_handler_3
      ParameterHandler prm;
      prm.enter_subsection ("Testing");
      prm.declare_entry ("string list",
                         "a",
                         Patterns::List(Patterns::Selection("a|b|c|d|e|f|g|h")),
                         "docs 1");
      prm.declare_entry ("int",
                         "1",
                         Patterns::Integer());
      prm.declare_entry ("double",
                         "3.1415926",
                         Patterns::Double(),
                         "docs 3");
      prm.leave_subsection ();
      
      prm.read_input("parameter_handler_3/prm");

                                       // now set some of the entries to
                                       // different values
      prm.enter_subsection ("Testing");
      prm.set ("string list", "a, c, b");
      prm.set ("int", "5");
      prm.set ("double", "2.71828");
      prm.leave_subsection ();
      
                                       // then write
      prm.print_parameters (logfile, ParameterHandler::Text);
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      
      return 1;
    }
  catch (...) 
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
  
  return 0;
}
