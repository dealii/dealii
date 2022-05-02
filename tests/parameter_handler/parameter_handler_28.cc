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
