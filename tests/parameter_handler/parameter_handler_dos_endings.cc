// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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


// ParameterHandler could not handle files with DOS line endings.

#include <deal.II/base/parameter_handler.h>

#include <sstream>
#include <string>

#include "../tests.h"


void
test()
{
  ParameterHandler foo;
  foo.enter_subsection("bar");
  foo.declare_entry("val", "1.0", dealii::Patterns::Double(), "");
  foo.leave_subsection();

  /*
   * At some point in the future deal.II may change its VCS settings so that
   * all files are automatically converted to Unix ('\n') line endings. Since
   * this bug involves files with DOS ('\r\n') line endings the parameter file
   * must be stored here as a string.
   */
  std::stringstream file_contents;
  file_contents << "# Note that this file has DOS (\\r\\n) line endings.\r\n";
  file_contents << "\r\n";
  file_contents << "subsection bar \r\n";
  file_contents << "  set val = 123.456\r\n";
  file_contents << "end\r\n";

  std::istringstream input_stream(file_contents.str());
  foo.parse_input(input_stream);

  foo.enter_subsection("bar");
  deallog << foo.get("val") << std::endl;
  foo.leave_subsection();
}

int
main()
{
  initlog();
  deallog.depth_console(0);

  test();

  return 0;
}
