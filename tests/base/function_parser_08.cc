// ---------------------------------------------------------------------
// $Id$
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



// test additional variables in function parser

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

  

  std::string function;
  std::map<std::string, double> constants;
  
  {
    function = "2*x + 0*y + depth";
    FunctionParser<2> fp;
    fp.initialize("x,y,notused,depth",
		  "2*x + depth", constants);
  
    Point<2> point(2.0, 3.0);
    std::vector<double> additional_vars(2);
    deallog << "Function " << "[" << function << "]" <<
          " @point " << "[" << point << "]" << " is " <<
      "[" <<  fp.value(point, 0, additional_vars) << "] with depth=0" << std::endl;

    additional_vars[1] = 3.0;
    
    deallog << "Function " << "[" << function << "]" <<
      " @point " << "[" << point << "]" << " is " <<
      "[" <<  fp.value(point, 0, additional_vars) << "] with depth=3" << std::endl;
  }

}
