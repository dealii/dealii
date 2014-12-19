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
// wrapper with respect to the ability to deal with units.
// because units are deprecated and muparser can not deal with this,
// recycle the test to use constants instead of units

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

  

  std::vector<std::string> function(1);
  std::map<std::string, double> constants;
  std::map<std::string, double> units;

  constants["PI"] = 3.141592654;
  constants["cm"] = 10;
  constants["m"] = 1000;

  Point<2> point(2.0, 3.0);

  //initialized with units
  FunctionParser<2> fp;
  function[0] = "x * cm + y * m + PI";
  fp.initialize(FunctionParser<2>::default_variable_names(),
                function, constants, units);

  deallog << "Function " << "[" << function[0] << "]" <<
          " @point " << "[" << point << "]" << " is " <<
          "[" <<  fp.value(point) << "]" << std::endl;

  //now initialize with a function
  //that's a string, not vector of
  //strings
  FunctionParser<2> fp4;
  fp4.initialize(FunctionParser<2>::default_variable_names(),
                 function[0], constants, units);

  deallog << "Function " << "[" << function[0] << "]" <<
          " @point " << "[" << point << "]" << " is " <<
          "[" <<  fp4.value(point) << "]" << std::endl;

  //now initialize a function without
  //units to check backwards
  //compatibility
  FunctionParser<2> fp2;
  function[0] = "x + y + PI";
  fp2.initialize(FunctionParser<2>::default_variable_names(),
                 function, constants);
  deallog << "Function " << "[" << function[0] << "]" <<
          " @point " << "[" << point << "]" << " is " <<
          "[" <<  fp2.value(point) << "]" << std::endl;



  //same as above but the function is
  //a string, not a vector
  FunctionParser<2> fp3;
  fp3.initialize(FunctionParser<2>::default_variable_names(),
                 function[0], constants);

  deallog << "Function " << "[" << function[0] << "]" <<
          " @point " << "[" << point << "]" << " is " <<
          "[" <<  fp3.value(point) << "]" << std::endl;
}




