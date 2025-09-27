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



// functionparser: check if you can change the expression in an existing object

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
  std::string variables = "s,t";

  FunctionParser<2> fp;
  fp.initialize("s,t", "s*t+1", constants);

  double value = fp.value(Point<2>(2.0, 2.5));
  Assert(std::abs(1.0 + 2.0 * 2.5 - value) < 1e-10, ExcMessage("wrong value"));

  std::vector<std::string> expressions;
  expressions.push_back("sin(2*mypi*x)+y");
  constants["mypi"] = numbers::PI;
  fp.initialize("x,y", expressions, constants);
  double value1 = fp.value(Point<2>(1.0, 2.5), 0);
  Assert(std::abs(2.5 - value1) < 1e-10, ExcMessage("wrong value"));
}
