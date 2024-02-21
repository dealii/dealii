// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check ParameterHandler::parse_input_from_xml. try to read a file with a
// parameter that does not satisfy its pattern

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

  // read from XML
  std::ifstream in(SOURCE_DIR "/prm/parameter_handler_read_xml_error_04.prm");
  try
    {
      prm.parse_input_from_xml(in);
    }
  catch (const ParameterHandler::ExcValueDoesNotMatchPattern &exc)
    {
      deallog << exc.get_exc_name() << std::endl;
      exc.print_info(deallog.get_file_stream());
    }

  return 0;
}
