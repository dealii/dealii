// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2018 by the deal.II authors
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

#include <deal.II/base/function_parser.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include <map>

#include "../tests.h"


int
main()
{
  initlog();

  // Define some constants that will be used by the function parser
  std::map<std::string, double> constants;
  constants["pi"] = numbers::PI;

  // Define the variables that will be used inside the expressions
  std::string variables = "x,y";

  // Define the expressions of the vector_valued function.
  std::vector<std::string> expressions;
  expressions.push_back("sin(2*pi*x)");
  expressions.push_back("cos(2*pi*y*x)");
  expressions.push_back("if(x<.5,y,exp(x))");
  expressions.push_back("if(x<.5,y,exp(x)*exp(x*y))");

  // Concatenate the declared expressions, to test the second way of
  // initializing
  std::string concatenated = "cos(pi*y)";
  // Now test each possibility
  for (unsigned int i = 1; i <= expressions.size(); ++i)
    {
      try
        {
          {
            FunctionParser<2> function(i);
            function.initialize(variables, expressions, constants);
            deallog << "Initialize Succeeded with dim = 2, " << i
                    << " components, " << expressions.size() << " expressions, "
                    << variables << " as variables." << std::endl;
          }
        }
      catch (...)
        {
          deallog << "Initialize Failed with dim = 2, " << i << " components, "
                  << expressions.size() << " expressions, " << variables
                  << " as variables." << std::endl;
        }

      try
        {
          {
            FunctionParser<2> function_bis(i);
            function_bis.initialize(variables, concatenated, constants);
            deallog << "Initialize Succeeded with dim = 2, " << i
                    << " components, " << concatenated << " as function and "
                    << variables << " as variables." << std::endl;
          }
        }
      catch (...)
        {
          deallog << "Initialize Failed with dim = 2, " << i << " components, "
                  << concatenated << " as function and " << variables
                  << " as variables." << std::endl;
        }

      concatenated += "; " + expressions[i - 1];
    }
}
