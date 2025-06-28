// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// This program tests the functionality of imposing a default expression
// when declaring a ParsedFunction. This functionality is tested for
// functions with 1, 2 and 3 components and when the default expression is
// not explicitly specified.

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include "deal.II/matrix_free/matrix_free.h"

#include "../tests.h"


void
test()
{
  // set up problem:
  ParameterHandler             prm;
  Functions::ParsedFunction<3> fct_1(1);
  Functions::ParsedFunction<3> fct_2(2);
  Functions::ParsedFunction<3> fct_3(3);

  std::string expression_1 = "1000";
  std::string expression_3 = "10;10;10";

  // Declare
  prm.enter_subsection("Function 1");
  fct_1.declare_parameters(prm, 1, expression_1);
  prm.leave_subsection();

  prm.enter_subsection("Function 2");
  fct_2.declare_parameters(prm, 2);
  prm.leave_subsection();

  prm.enter_subsection("Function 3");
  fct_3.declare_parameters(prm, 3, expression_3);
  prm.leave_subsection();


  // Parse
  prm.enter_subsection("Function 1");
  fct_1.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Function 2");
  fct_2.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Function 3");
  fct_3.parse_parameters(prm);
  prm.leave_subsection();


  // Dummy point
  const Point<3> p;

  const double   result_1 = fct_1.value(p);
  Vector<double> result_2(2);
  Vector<double> result_3(3);

  fct_2.vector_value(p, result_2);
  fct_3.vector_value(p, result_3);

  deallog << "Function 1 '" << expression_1 << "' is " << result_1 << std::endl;

  deallog << "Function 2 '' is " << result_2 << std::endl;

  deallog << "Function 3 '" << expression_3 << "' is " << result_3 << std::endl;
}

int
main()
{
  initlog();

  test();
}
