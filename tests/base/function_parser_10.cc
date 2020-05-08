// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// test rand function and expression only constructor

#include <deal.II/base/function_parser.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include <map>

#include "../tests.h"


double
eval(const std::string &exp)
{
  FunctionParser<2> fp(exp);

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
