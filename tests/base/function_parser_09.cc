// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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


bool
satisfies_randomness(std::vector<double> &v)
{
  // check that consecutive numbers are not the same
  for (unsigned int i = 1; i < v.size(); ++i)
    if (v[i - 1] == v[i])
      return false;
  // Check that the numbers are between 0 and 1
  return std::all_of(v.begin(), v.end(), [](double n) {
    return (n >= 0.) && (n <= 1.);
  });
}


int
main()
{
  initlog();

  std::vector<double> rands{eval("rand_seed(10)"),
                            eval("rand()"),
                            eval("rand_seed(10)"),
                            eval("rand_seed(10)")};

  if (satisfies_randomness(rands))
    deallog << "OK" << std::endl;
}
