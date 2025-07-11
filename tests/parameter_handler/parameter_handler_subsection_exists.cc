// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test ParameterHandler::subsection_path_exists().

#include <deal.II/base/parameter_handler.h>

#include <iomanip>

#include "../tests.h"

int
main()
{
  initlog();

  ParameterHandler prm;
  prm.enter_subsection("Section 1");
  {
    prm.declare_entry("string list",
                      "a",
                      Patterns::List(Patterns::Selection("a|b|c|d|e|f|g|h")),
                      "docs 1");
    prm.declare_entry("int", "1", Patterns::Integer());
    prm.declare_entry("double", "3.1415926", Patterns::Double(), "docs 3");
  }
  prm.leave_subsection();

  prm.enter_subsection("Section 2");
  {
    prm.declare_entry("int", "1", Patterns::Integer());
    prm.declare_entry("double", "3.1415926", Patterns::Double(), "docs 3");

    prm.enter_subsection("string list");
    {
      prm.declare_entry("string list",
                        "a",
                        Patterns::List(Patterns::Selection("a|b|c|d|e|f|g|h")),
                        "docs 1");
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();

  prm.print_parameters(deallog.get_file_stream(),
                       ParameterHandler::PRM |
                         ParameterHandler::KeepDeclarationOrder);

  deallog << std::boolalpha;
  deallog << "Subsection \"Section 3\" of root exists: "
          << prm.subsection_path_exists({"Section 3"}) << std::endl;
  deallog << "Subsection \"Section 1.string list\" of root exists: "
          << prm.subsection_path_exists({"Section 1", "string list"})
          << std::endl;

  deallog << "Subsection \"Section 2.string list\" of root exists: "
          << prm.subsection_path_exists({"Section 2", "string list"})
          << std::endl;

  prm.enter_subsection("Section 2");
  {
    deallog << "Subsection \"string list\" of \"Section 2\" exists: "
            << prm.subsection_path_exists({"string list"}) << std::endl;
  }
  prm.leave_subsection();

  deallog
    << "Subsection \"Section 2.string list\" still exists after coming back to root: "
    << prm.subsection_path_exists({"Section 2", "string list"}) << std::endl;

  return 0;
}
