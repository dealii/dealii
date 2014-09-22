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



// functionparser: check if you can change the expression in an existing object

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
  initlog();

  // Define some constants that will be used by the function parser
  std::map<std::string, double> constants;
  constants["pi"] = numbers::PI;

  // Define the variables that will be used inside the expressions
  std::string variables = "s,t";

  FunctionParser<2> fp;
  fp.initialize("s,t", "s*t+1", constants);

  double value = fp.value(Point<2>(2.0,2.5));
  Assert(abs(1.0+2.0*2.5 - value) < 1e-10, ExcMessage("wrong value"));
  
  std::vector<std::string> expressions;
  expressions.push_back("sin(2*mypi*x)+y");
  constants["mypi"] = numbers::PI;
  fp.initialize("x,y", expressions, constants);
  double value1 = fp.value(Point<2>(1.0,2.5), 0);
  Assert(abs(2.5 - value1) < 1e-10, ExcMessage("wrong value"));
  
}




