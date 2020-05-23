// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2019 by the deal.II authors
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



// This program tests the functionality of the function parser
// wrapper.

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_function_parser.h>

#include "../tests.h"


// Test initializtaion and evaluation of timedependent function in 2d.
int
main()
{
  initlog();

  // Define some constants that will be used by the function parser
  std::map<std::string, double> constants;
  constants["pi"] = numbers::PI;

  // Define the variables that will be used inside the expressions
  std::string variables = "x,y,t";

  // Define the expressions of the vector_valued function.
  std::vector<std::string> expressions;
  expressions.push_back("cos(2*pi*(x*y+t))");
  expressions.push_back("sin(2*pi*(x*y+t))");
  expressions.push_back("-sin(2*pi*(x*y+t))");
  expressions.push_back("cos(2*pi*(x*y+t))");

  // Test time dependent function
  try
    {
      {
        TensorFunctionParser<2, 2> tensor_function;
        tensor_function.initialize(variables,
                                   expressions,
                                   constants,
                                   /* time dependent */ true);

        deallog << "Initialize Succeeded with dim = 2, rank = 2, "
                << expressions.size() << " expressions, " << variables
                << " as variables." << std::endl;
      }
    }
  catch (...)
    {
      deallog << "Initialization or Evaluation Failed with dim = 2, rank = 2, "
              << expressions.size() << " expressions, " << variables
              << " as variables." << std::endl;
    }
}
