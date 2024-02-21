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



// test ParameterHandler::print_parameters for different formats in
// combination with Short

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

int
main()
{
  initlog();

  double a = 1.0;
  double b = 2.0;

  ParameterHandler prm;
  prm.add_parameter("A double", a, "Documentation for a.");
  prm.enter_subsection("Section two");
  prm.add_parameter("Another double", b, "Documentation for b.");
  prm.leave_subsection();

  const auto style = ParameterHandler::Short;

  deallog << "ShortPRM format: " << std::endl
          << "========================================" << std::endl;
  prm.print_parameters("output.prm", style);
  cat_file("output.prm");

  deallog << "ShortXML format: " << std::endl
          << "========================================" << std::endl;
  prm.print_parameters("output.xml", style);
  cat_file("output.xml");

  deallog << "ShortJSON format: " << std::endl
          << "========================================" << std::endl;
  prm.print_parameters("output.json", style);
  cat_file("output.json");

  deallog << "ShortLaTeX format: " << std::endl
          << "========================================" << std::endl;
  prm.print_parameters("output.tex", style);
  cat_file("output.tex");

  deallog << "Short | Description format: " << std::endl
          << "========================================" << std::endl;
  prm.print_parameters("output.dsc", ParameterHandler::Description | style);
  cat_file("output.dsc");
}
