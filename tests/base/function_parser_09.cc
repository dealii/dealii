// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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


// test rand function

#include <deal.II/base/function_parser.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include <map>

#include "../tests.h"


double
eval(const std::string &exp)
{
  std::string                   variables = "x,y";
  std::map<std::string, double> constants;

  FunctionParser<2> fp(1);
  fp.initialize(variables, exp, constants);

  Point<2> p;

  return fp.value(p);
}


int
main()
{
  initlog();

  double random = eval("rand()"); // random seed

  if (0.0 <= random && random <= 1.0)
    deallog << "OK" << std::endl;

  deallog << eval("rand_seed(10)") << std::endl;
}
