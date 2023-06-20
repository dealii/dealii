// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2023 by the deal.II authors
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
