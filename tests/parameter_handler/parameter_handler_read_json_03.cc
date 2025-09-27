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



// check ParameterHandler::parse_input() for json file

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
  int int5 = 0;

  ParameterHandler prm;
  prm.add_parameter("int1", int1);
  prm.add_parameter("int2", int2);

  prm.enter_subsection("dummy");
  prm.add_parameter("int3", int3);
  prm.add_parameter("int4", int4);
  prm.leave_subsection();

  prm.enter_subsection("test_underscore");
  prm.add_parameter("int_5", int5);
  prm.leave_subsection();

  std::string source   = SOURCE_DIR;
  std::string filename = source + "/prm/parameter_handler_read_json_03.json";
  prm.parse_input(filename, "", true, true);

  AssertDimension(int1, 1);
  AssertDimension(int2, 2);
  AssertDimension(int3, 3);
  AssertDimension(int4, 4);
  AssertDimension(int5, 10);

  return 0;
}
