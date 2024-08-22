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



// Like parameter_handler_read_json_05 but output with
// ParameterHandler::OutputStyle::KeepOnlyChanged. This is achieved
// by adding two parameters with default values, altering one of
// the value of the parameters and printing the parameters
// using the ParameterHandler::OutputStyle::KeepOnlyChanged.

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


int
main()
{
  initlog();

  ParameterHandler prm;

  double test_0 = 0;
  double test_1 = 1;

  // test if underscore can be parsed
  prm.add_parameter("test 0", test_0);
  prm.add_parameter("test 1", test_1);

  std::string source   = SOURCE_DIR;
  std::string filename = source + "/prm/parameter_handler_read_json_04.json";

  std::ifstream file;
  file.open(filename);
  prm.parse_input_from_json(file, true);

  prm.print_parameters(deallog.get_file_stream(),
                       ParameterHandler::OutputStyle::ShortJSON |
                         ParameterHandler::OutputStyle::KeepOnlyChanged);

  return 0;
}
