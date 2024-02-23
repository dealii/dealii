// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test ParameterHandler::print_parameters for different formats

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

int
main()
{
  initlog();

  double      a         = 1.0;
  bool        some_bool = true;
  std::string a_string  = "Ciao";

  ParameterHandler prm;
  prm.enter_subsection("Testing");
  prm.add_parameter("A double", a);
  prm.add_parameter("A bool", some_bool);
  prm.add_parameter("A string", a_string);
  prm.leave_subsection();

  deallog << "PRM format: " << std::endl
          << "========================================" << std::endl;
  prm.print_parameters("output.prm");
  cat_file("output.prm");

  deallog << "XML format: " << std::endl
          << "========================================" << std::endl;
  prm.print_parameters("output.xml");
  cat_file("output.xml");

  deallog << "JSON format: " << std::endl
          << "========================================" << std::endl;
  prm.print_parameters("output.json");
  cat_file("output.json");

  deallog << "LaTeX format: " << std::endl
          << "========================================" << std::endl;
  prm.print_parameters("output.tex");
  cat_file("output.tex");

  deallog << "Description format: " << std::endl
          << "========================================" << std::endl;
  prm.print_parameters("output.dsc", ParameterHandler::Description);
  cat_file("output.dsc");
}
