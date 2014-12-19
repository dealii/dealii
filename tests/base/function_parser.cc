// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// This program tests the functionality of the function parser
// wrapper.

#include "../tests.h"
#include <fstream>
#include <iomanip>
#include <map>
#include <deal.II/base/logstream.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/function_parser.h>


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

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
  std::string concatenated="cos(pi*y)";
  // Now test each possibility
  for (unsigned int i=1; i<=expressions.size(); ++i)
    {
      try
        {
          {
            FunctionParser<2> function(i);
            function.initialize(variables,
                                expressions,
                                constants);
            deallog << "Initialize Succeded with dim = 2, "
                    << i << " components, "
                    << expressions.size() << " expressions, "
                    << variables << " as variables." << std::endl;
          }
        }
      catch (...)
        {
          deallog << "Initialize Failed with dim = 2, "
                  << i << " components, "
                  << expressions.size() << " expressions, "
                  << variables << " as variables." << std::endl;
        }

      try
        {
          {
            FunctionParser<2> function_bis(i);
            function_bis.initialize(variables,
                                    concatenated,
                                    constants);
            deallog << "Initialize Succeded with dim = 2, "
                    << i << " components, "
                    << concatenated << " as function and "
                    << variables << " as variables." << std::endl;
          }
        }
      catch (...)
        {
          deallog << "Initialize Failed with dim = 2, "
                  << i << " components, "
                  << concatenated  << " as function and "
                  << variables << " as variables." << std::endl;
        }

      concatenated += "; " + expressions[i-1];
    }
}




