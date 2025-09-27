// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Document a bug that rand_seed() in the FunctionParser always
// returns the same value.

#include <deal.II/base/function_parser.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include <iostream>
#include <map>
#include <random>

#include "../tests.h"


// This is copied from the FunctionParser implementation
double
mu_rand_seed(double seed)
{
  std::uniform_real_distribution<> uniform_distribution(0., 1.);

  static std::map<double, std::mt19937> rng_map;

  if (rng_map.find(seed) == rng_map.end())
    rng_map[seed] = std::mt19937(static_cast<unsigned int>(seed));

  return uniform_distribution(rng_map[seed]);
}


int
main()
{
  initlog();
  double                        seed      = 1;
  std::string                   variables = "x,y";
  std::map<std::string, double> constants;

  FunctionParser<2> fp(1);
  fp.initialize(variables, "rand_seed(1)", constants);
  Point<2> p;

  for (int i = 0; i < 10; ++i)
    {
      const double manual  = mu_rand_seed(seed);
      const double from_fp = fp.value(p);
      AssertThrow(abs(manual - from_fp) < 1e-16, ExcInternalError());
    }

  deallog << "OK" << std::endl;
}
