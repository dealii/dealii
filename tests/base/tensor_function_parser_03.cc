// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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

  Point<2> p(0.5, 0.5); // evaluation point

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

        for (unsigned int t = 0; t < 5; ++t)
          {
            tensor_function.set_time(t);
            deallog << "Value at p = (0.5,0.5):   " << std::endl
                    << tensor_function.value(p)[0][0] << "   "
                    << tensor_function.value(p)[0][1] << std::endl
                    << tensor_function.value(p)[1][0] << "   "
                    << tensor_function.value(p)[1][1] << "   at time t = " << t
                    << std::endl;
          }
      }
    }
  catch (...)
    {
      deallog << "Initialization or Evaluation Failed with dim = 2, rank = 2, "
              << expressions.size() << " expressions, " << variables
              << " as variables." << std::endl;
    }
}
