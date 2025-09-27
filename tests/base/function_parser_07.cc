// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// like _06, but ensure that we can not only parse "if(cond,yes,no)" but also
// "if (cond,yes,no)", i.e., including the whitespace between function and
// argument list

#include <deal.II/base/function_parser.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include <map>

#include "../tests.h"


void
eval(const std::string &exp, const Point<2> &p, double expected)
{
  std::string                   variables = "x,y";
  std::map<std::string, double> constants;

  FunctionParser<2> fp(1);
  fp.initialize(variables, exp, constants);

  double result = fp.value(p);
  deallog << "'" << exp << "' @ " << p << " is " << result << " ( expected "
          << expected << " )" << std::endl;
  if (fabs(result - expected) > 1e-10)
    deallog << "ERROR!" << std::endl;
}


void
test()
{
  eval("if   (x<0.0,0.0,if(x>1.0,1.0,x))", Point<2>(0.5, 0.0), 0.5);
  eval("if   (x<0.0,0.0,if(x>1.0,1.0,x))", Point<2>(-2.0, 0.0), 0.0);
  eval("if   (x<0.0,0.0,if(x>1.0,1.0,x))", Point<2>(42.0, 0.0), 1.0);

  eval("if   (x<1.0 | y<1.0,0,y)", Point<2>(0.5, 1.5), 0.0);
  eval("if   (x<1.0 | y<1.0,0,y)", Point<2>(1.5, 1.5), 1.5);
  eval("if   (x<1.0 | y<1.0,0,y)", Point<2>(1.5, 0.5), 0.0);
  eval("if   (x<1.0 | y<1.0,0,y)", Point<2>(0.5, -2), 0.0);

  eval("if   (x<1.0 & y<1.0,0,y)", Point<2>(1.5, -2.0), -2.0);

  double x, y;
  x = 1.0;
  y = -3.1;
  eval("atan2 (x,y)", Point<2>(x, y), atan2(x, y));
  x = -1.0;
  y = 3.1;
  eval("atan2 (x,y)", Point<2>(x, y), atan2(x, y));

  eval("if   (x==1.0,0,y)", Point<2>(1.0, -2.0), 0.0);
  eval("if   (x==1.0,0,y)", Point<2>(1.1, -2.0), -2.0);

  eval("int (2.1)", Point<2>(1.1, -2.0), 2.0);
  eval("int (-3.8)", Point<2>(1.1, -2.0), -4.0);

  eval("abs (-2.3)", Point<2>(0, 0), 2.3);
  eval("acos (0.5)", Point<2>(0, 0), acos(0.5));
  eval("acosh (1.5)", Point<2>(0, 0), acosh(1.5));
  eval("asin (0.5)", Point<2>(0, 0), asin(0.5));
  eval("asinh (0.5)", Point<2>(0, 0), asinh(0.5));
  eval("atan (0.5)", Point<2>(0, 0), atan(0.5));
  eval("atanh (0.5)", Point<2>(0, 0), atanh(0.5));
  eval("ceil (0.5)", Point<2>(0, 0), ceil(0.5));
  eval("cos (0.5)", Point<2>(0, 0), cos(0.5));
  eval("cosh (0.5)", Point<2>(0, 0), cosh(0.5));
  eval("cot (0.5)", Point<2>(0, 0), 1.0 / tan(0.5));
  eval("csc (0.5)", Point<2>(0, 0), 1.0 / sin(0.5));
  eval("exp (0.5)", Point<2>(0, 0), exp(0.5));
  eval("floor (0.5)", Point<2>(0, 0), floor(0.5));
  eval("log (0.5)", Point<2>(0, 0), log(0.5));
  eval("log10 (0.5)", Point<2>(0, 0), log(0.5) / log(10.0));
  eval("sec (0.5)", Point<2>(0, 0), 1.0 / cos(0.5));
  eval("sin (0.5)", Point<2>(0, 0), sin(0.5));
  eval("sinh (0.5)", Point<2>(0, 0), sinh(0.5));
  eval("sqrt (0.5)", Point<2>(0, 0), sqrt(0.5));
  eval("tan (0.5)", Point<2>(0, 0), tan(0.5));
  eval("tanh (0.5)", Point<2>(0, 0), tanh(0.5));
}

int
main()
{
  initlog();

  test();
}
