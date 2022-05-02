// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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
                       ParameterHandler::Text |
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
