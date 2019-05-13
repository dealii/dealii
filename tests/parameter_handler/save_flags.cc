// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



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
