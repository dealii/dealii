//----------------------------  parameter_handler_5.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  parameter_handler_5.cc  ---------------------------


// ParameterHandler::declare_entry did not allow to redeclare an
// entry. make sure this works now

#include <base/logstream.h>
#include <base/parameter_handler.h>
#include <fstream>
#include <iostream>


int main () 
{
  try
    {
      std::ofstream logfile("parameter_handler_5.output");
      deallog.attach(logfile);
      deallog.depth_console(0);

      ParameterHandler prm;
      prm.declare_entry ("int",
                         "1",
                         Patterns::Integer());
      prm.declare_entry ("int",
                         "2",
                         Patterns::Integer());

      prm.print_parameters (logfile, ParameterHandler::Text);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      
      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
  
  return 0;
}
