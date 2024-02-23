// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that ParameterHandler::print_parameters() saves and resets
// the iostream flags of the stream it writes to

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


int
main()
{
  initlog();

  ParameterHandler prm;
  prm.declare_entry("int1", "1", Patterns::Integer(), "doc 1");
  prm.declare_entry("int2", "2", Patterns::Integer(), "doc 2");
  prm.enter_subsection("ss1");
  {
    prm.declare_entry("double 1", "1.234", Patterns::Double(), "doc 3");

    prm.enter_subsection("ss2");
    {
      prm.declare_entry("double 2", "4.321", Patterns::Double(), "doc 4");
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

  // set a special fill char and verify that it is being used
  deallog.get_file_stream().fill('x');
  deallog.get_file_stream().width(15);
  deallog.get_file_stream() << std::left << 42 << std::endl;

  // now let ParameterHandler output its state
  prm.print_parameters(deallog.get_file_stream(),
                       ParameterHandler::Description);

  // verify that the special fill char is still available (i.e., that
  // print_parameters() has saved and restored the stream flags)
  deallog.get_file_stream().width(15);
  deallog.get_file_stream() << std::left << 42 << std::endl;

  return 0;
}
