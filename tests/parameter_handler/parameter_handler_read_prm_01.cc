// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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



// check ParameterHandler::parse_input for prm file

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


int
main()
{
  initlog();

  // default values
  int int1 = 0;
  int int2 = 0;
  int int3 = 0;
  int int4 = 0;

  ParameterHandler prm;
  prm.add_parameter("int1", int1);
  prm.add_parameter("int2", int2);

  prm.enter_subsection("dummy");
  prm.add_parameter("int3", int3);
  prm.add_parameter("int4", int4);
  prm.leave_subsection();

  // read from json
  std::string source   = SOURCE_DIR;
  std::string filename = source + "/prm/parameter_handler_read_prm_01.prm";
  prm.parse_input(filename, "", true, true);

  AssertDimension(int1, 1);
  AssertDimension(int2, 2);
  AssertDimension(int3, 3);
  AssertDimension(int4, 4);

  return 0;
}
