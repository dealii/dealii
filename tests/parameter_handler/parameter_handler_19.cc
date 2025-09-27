// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// ParameterHandler does not complain if you parse input that doesn't close
// all subsections.


#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

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
  catch (const std::exception &e)
    {
      deallog << "Exception " << e.what() << std::endl;
    }



  deallog << std::endl << "* parse_input with missing 'end':" << std::endl;

  try
    {
      std::string s = "set dim=2\nsubsection test\n\n"; // note: missing "end"
      prm.parse_input_from_string(s.c_str());
    }
  catch (const ParameterHandler::ExcUnbalancedSubsections &exc)
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
  catch (const std::exception &e)
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
    catch (const ParameterHandler::ExcUnbalancedSubsections &exc)
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
