// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2017 by the deal.II authors
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

// ParameterHandler does not complain if you parse input that doesn't close
// all subsections.

#include "../tests.h"
#include <deal.II/base/parameter_handler.h>

void
check()
{
  ParameterHandler prm;
  prm.declare_entry("dim", "3", Patterns::Integer());
  prm.enter_subsection("test");
  prm.declare_entry("x", "1", Patterns::Integer());
  prm.leave_subsection();
  prm.enter_subsection("test2");
  prm.declare_entry("y", "1", Patterns::Integer());
  prm.leave_subsection();

  deallog << "* no subsection to leave: " << std::endl;
  try
    {
      prm.leave_subsection();
    }
  catch(std::exception& e)
    {
      deallog << "Exception " << e.what() << std::endl;
    }

  deallog << std::endl << "* parse_input with missing 'end':" << std::endl;

  try
    {
      std::string s = "set dim=2\nsubsection test\n\n"; // note: missing "end"
      prm.parse_input_from_string(s.c_str());
    }
  catch(const ParameterHandler::ExcUnbalancedSubsections& exc)
    {
      deallog << exc.get_exc_name() << std::endl;
      exc.print_info(deallog.get_file_stream());
    }
  deallog << std::endl;

  // make sure parse_input resets the current path:
  try
    {
      prm.leave_subsection();
      deallog << "error, why could we leave a subsection?" << std::endl;
    }
  catch(std::exception& e)
    {
      deallog << "Exception " << e.what() << std::endl;
    }

  deallog << std::endl
          << "* Check non empty path before parse_input()" << std::endl;

  {
    prm.enter_subsection("test");
    std::string s = "set x=5\n";
    prm.parse_input_from_string(s.c_str());
    prm.leave_subsection();
  }

  deallog << std::endl
          << "* Check parse_input() catches messing with path:" << std::endl;
  {
    prm.enter_subsection("test");
    try
      {
        std::string s = "end\nsubsection test2\nset y=7\n";
        prm.parse_input_from_string(s.c_str());
      }
    catch(const ParameterHandler::ExcUnbalancedSubsections& exc)
      {
        deallog << exc.get_exc_name() << std::endl;
        exc.print_info(deallog.get_file_stream());
      }
    prm.leave_subsection();
  }
}

int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  check();

  return 0;
}
