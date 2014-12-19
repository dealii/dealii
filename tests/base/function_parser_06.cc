// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2014 by the deal.II authors
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

// various checks (compatibility functionparser -> muparser)


#include "../tests.h"
#include <fstream>
#include <iomanip>
#include <map>
#include <deal.II/base/logstream.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/function_parser.h>


void eval(const std::string & exp, const Point<2> & p, double expected)
{
  std::string variables = "x,y";
  std::map<std::string,double> constants;

  FunctionParser<2> fp(1);
  fp.initialize(variables,
		exp,
		constants);

  double result = fp.value(p);
  deallog << "'" << exp << "' @ " << p << " is " << result
	  << " ( expected " << expected << " )" << std::endl;
  if (fabs(result-expected)>1e-10)
    deallog << "ERROR!" << std::endl;
  
  
  
}


void test()
{
  eval("if(x<0.0,0.0,if(x>1.0,1.0,x))",Point<2>(0.5,0.0), 0.5);
  eval("if(x<0.0,0.0,if(x>1.0,1.0,x))",Point<2>(-2.0,0.0), 0.0);
  eval("if(x<0.0,0.0,if(x>1.0,1.0,x))",Point<2>(42.0,0.0), 1.0);

  eval("if(x<1.0 | y<1.0,0,y)",Point<2>(0.5,1.5), 0.0);
  eval("if(x<1.0 | y<1.0,0,y)",Point<2>(1.5,1.5), 1.5);
  eval("if(x<1.0 | y<1.0,0,y)",Point<2>(1.5,0.5), 0.0);
  eval("if(x<1.0 | y<1.0,0,y)",Point<2>(0.5,-2), 0.0);
  
  eval("if(x<1.0 & y<1.0,0,y)",Point<2>(1.5,-2.0), -2.0);

  double x,y;
  x=1.0;
  y=-3.1;
  eval("atan2(x,y)",Point<2>(x,y), atan2(x,y));
  x=-1.0;
  y=3.1;
  eval("atan2(x,y)",Point<2>(x,y), atan2(x,y));

  eval("if(x==1.0,0,y)",Point<2>(1.0,-2.0), 0.0);  
  eval("if(x==1.0,0,y)",Point<2>(1.1,-2.0), -2.0);

  eval("int(2.1)",Point<2>(1.1,-2.0), 2.0);
  eval("int(-3.8)",Point<2>(1.1,-2.0), -4.0);

  eval("abs(-2.3)",Point<2>(0,0), 2.3);
  eval("acos(0.5)",Point<2>(0,0), acos(0.5));
  eval("acosh(0.5)",Point<2>(0,0), acosh(0.5));
  eval("asin(0.5)",Point<2>(0,0), asin(0.5));
  eval("asinh(0.5)",Point<2>(0,0), asinh(0.5));
  eval("atan(0.5)",Point<2>(0,0), atan(0.5));
  eval("atanh(0.5)",Point<2>(0,0), atanh(0.5));
  eval("ceil(0.5)",Point<2>(0,0), ceil(0.5));
  eval("cos(0.5)",Point<2>(0,0), cos(0.5));
  eval("cosh(0.5)",Point<2>(0,0), cosh(0.5));
  eval("cot(0.5)",Point<2>(0,0), 1.0/tan(0.5));
  eval("csc(0.5)",Point<2>(0,0), 1.0/sin(0.5));
  eval("exp(0.5)",Point<2>(0,0), exp(0.5));
  eval("floor(0.5)",Point<2>(0,0), floor(0.5));
  eval("log(0.5)",Point<2>(0,0), log(0.5));
  eval("log10(0.5)",Point<2>(0,0), log(0.5)/log(10.0));
  eval("sec(0.5)",Point<2>(0,0), 1.0/cos(0.5));
  eval("sin(0.5)",Point<2>(0,0), sin(0.5));
  eval("sinh(0.5)",Point<2>(0,0), sinh(0.5));
  eval("sqrt(0.5)",Point<2>(0,0), sqrt(0.5));
  eval("tan(0.5)",Point<2>(0,0), tan(0.5));
  eval("tanh(0.5)",Point<2>(0,0), tanh(0.5));


}

int main()
{
  initlog();

  test();
}
