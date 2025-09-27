// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check ParameterHandler::print_parameters (..., JSON). have a few
// names that contain all sorts of weird (for JSON) characters and
// aliased parameters.

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


int
main()
{
  initlog();

  ParameterHandler prm;
  prm.declare_entry("int1", "1", Patterns::Integer(), "doc 1");
  prm.declare_entry("int2", "2", Patterns::Integer(), "doc 2");
  prm.declare_alias("int2", "int2_alias");
  prm.enter_subsection("ss1");
  {
    prm.declare_entry("double 1", "1.234", Patterns::Double(), "doc 3");

    prm.enter_subsection("ss2");
    {
      prm.declare_entry("double 2", "4.321", Patterns::Double(), "doc 4");
      prm.declare_alias("double 2", "double 2 alias");
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();

  // things with strange characters
  prm.enter_subsection("Testing%testing");
  {
    prm.declare_entry("string&list",
                      "< & > ; /",
                      Patterns::Anything(),
                      "docs 1");
    prm.declare_entry("int*int", "2", Patterns::Integer());
    prm.declare_entry("double+double",
                      "6.1415926",
                      Patterns::Double(),
                      "docs 3");
  }
  prm.leave_subsection();


  prm.print_parameters(deallog.get_file_stream(), ParameterHandler::JSON);
  deallog.get_file_stream() << std::endl;

  prm.print_parameters(deallog.get_file_stream(),
                       ParameterHandler::JSON | ParameterHandler::Short);
  deallog.get_file_stream() << std::endl;


  return 0;
}
