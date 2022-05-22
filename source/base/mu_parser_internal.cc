// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2022 by the deal.II authors
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

#include <deal.II/base/mu_parser_internal.h>
#include <deal.II/base/thread_management.h>

#include <cmath>
#include <ctime>
#include <map>
#include <mutex>
#include <random>
#include <vector>


DEAL_II_NAMESPACE_OPEN



#ifdef DEAL_II_WITH_MUPARSER

namespace internal
{
  namespace FunctionParser
  {
    int
    mu_round(double val)
    {
      return static_cast<int>(val + ((val >= 0.0) ? 0.5 : -0.5));
    }

    double
    mu_if(double condition, double thenvalue, double elsevalue)
    {
      if (mu_round(condition) != 0)
        return thenvalue;
      else
        return elsevalue;
    }

    double
    mu_or(double left, double right)
    {
      return static_cast<double>((mu_round(left) != 0) ||
                                 (mu_round(right) != 0));
    }

    double
    mu_and(double left, double right)
    {
      return static_cast<double>((mu_round(left) != 0) &&
                                 (mu_round(right) != 0));
    }

    double
    mu_int(double value)
    {
      return static_cast<double>(mu_round(value));
    }

    double
    mu_ceil(double value)
    {
      return std::ceil(value);
    }

    double
    mu_floor(double value)
    {
      return std::floor(value);
    }

    double
    mu_cot(double value)
    {
      return 1.0 / std::tan(value);
    }

    double
    mu_csc(double value)
    {
      return 1.0 / std::sin(value);
    }

    double
    mu_sec(double value)
    {
      return 1.0 / std::cos(value);
    }

    double
    mu_log(double value)
    {
      return std::log(value);
    }

    double
    mu_pow(double a, double b)
    {
      return std::pow(a, b);
    }

    double
    mu_erf(double value)
    {
      return std::erf(value);
    }

    double
    mu_erfc(double value)
    {
      return std::erfc(value);
    }

    // returns a random value in the range [0,1] initializing the generator
    // with the given seed
    double
    mu_rand_seed(double seed)
    {
      static std::mutex           rand_mutex;
      std::lock_guard<std::mutex> lock(rand_mutex);

      std::uniform_real_distribution<> uniform_distribution(0., 1.);

      // for each seed a unique random number generator is created,
      // which is initialized with the seed itself
      static std::map<double, std::mt19937> rng_map;

      if (rng_map.find(seed) == rng_map.end())
        rng_map[seed] = std::mt19937(static_cast<unsigned int>(seed));

      return uniform_distribution(rng_map[seed]);
    }

    // returns a random value in the range [0,1]
    double
    mu_rand()
    {
      static std::mutex                rand_mutex;
      std::lock_guard<std::mutex>      lock(rand_mutex);
      std::uniform_real_distribution<> uniform_distribution(0., 1.);
      const unsigned int  seed = static_cast<unsigned long>(std::time(nullptr));
      static std::mt19937 rng(seed);
      return uniform_distribution(rng);
    }

    std::vector<std::string> function_names = {
      // functions predefined by muparser
      "sin",
      "cos",
      "tan",
      "asin",
      "acos",
      "atan",
      "sinh",
      "cosh",
      "tanh",
      "asinh",
      "acosh",
      "atanh",
      "atan2",
      "log2",
      "log10",
      "log",
      "ln",
      "exp",
      "sqrt",
      "sign",
      "rint",
      "abs",
      "min",
      "max",
      "sum",
      "avg",
      // functions we define ourselves above
      "if",
      "int",
      "ceil",
      "cot",
      "csc",
      "floor",
      "sec",
      "pow",
      "erf",
      "erfc",
      "rand",
      "rand_seed"};

  } // namespace FunctionParser

} // namespace internal
#endif



DEAL_II_NAMESPACE_CLOSE
