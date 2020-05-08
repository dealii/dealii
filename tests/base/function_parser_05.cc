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


#include <deal.II/base/function_parser.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include <map>

#include "../tests.h"

// simple function parser example (this is the basis for the example code in the
// documentation)

void
test1()
{
  // set up problem:
  std::string                   variables  = "x,y";
  std::string                   expression = "cos(x)+sqrt(y)";
  std::map<std::string, double> constants;

  // FunctionParser with 2 variables and 1 component:
  FunctionParser<2> fp(1);
  fp.initialize(variables, expression, constants);

  // Point at which we want to evaluate the function
  Point<2> point(0.0, 4.0);

  // evaluate the expression at 'point':
  double result = fp.value(point);

  deallog << "Function '" << expression << "'"
          << " @ " << point << " is " << result << std::endl;
}


void
test2()
{
  // Define some constants that will be used by the function parser
  std::map<std::string, double> constants;
  constants["pi"] = numbers::PI;

  // Define the variables that will be used inside the expressions
  std::string variables = "x,y,z";

  // Define the expressions of the individual components of a
  // vector valued function with two components:
  std::vector<std::string> expressions(2);
  expressions[0] = "sin(2*pi*x)+sinh(pi*z)";
  expressions[1] = "sin(2*pi*y)*exp(x^2)";

  // function parser with 3 variables and 2 components
  FunctionParser<3> vector_function(2);

  // And populate it with the newly created objects.
  vector_function.initialize(variables, expressions, constants);

  // Point at which we want to evaluate the function
  Point<3> point(0.0, 1.0, 1.0);

  // This Vector will store the result
  Vector<double> result(2);

  // Fill 'result' by evaluating the function
  vector_function.vector_value(point, result);

  // We can also only evaluate the 2nd component:
  double c = vector_function.value(point, 1);

  // Output the evaluated function
  deallog << "Function '" << expressions[0] << "," << expressions[1] << "'"
          << " @ " << point << " is " << std::flush;
  result.print(deallog.get_file_stream(), 4, false);
  deallog.get_file_stream() << "DEAL::" << std::endl;
}


int
main()
{
  initlog();

  test1();
  test2();
}
