// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2017 by the deal.II authors
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



// test ParameterHandler::set(., bool) which was broken, see bug #49

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


int
main()
{
  try
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);

      // same as parameter_handler_3
      ParameterHandler prm;
      prm.enter_subsection("Testing");
      prm.declare_entry("bool", "true", Patterns::Bool(), "docs 1");
      prm.leave_subsection();

      prm.parse_input(SOURCE_DIR "/prm/parameter_handler_6_bool.prm");

      // now set the parameter to a different
      // value
      prm.enter_subsection("Testing");
      prm.set("bool", true);
      prm.leave_subsection();

      // then write
      prm.print_parameters(logfile, ParameterHandler::Text);

      // and do it again with the opposite
      // value
      prm.enter_subsection("Testing");
      prm.set("bool", false);
      prm.leave_subsection();

      // then write
      prm.print_parameters(logfile, ParameterHandler::Text);
    }
  catch (std::exception &exc)
    {
      deallog << std::endl
              << std::endl
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
      deallog << std::endl
              << std::endl
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
