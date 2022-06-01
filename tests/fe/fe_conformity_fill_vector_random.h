// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#pragma once

#include <deal.II/lac/vector.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <typeinfo>
#include <vector>

namespace FEConforimityTest
{
  using namespace dealii;

  /*
   * Generate a random double number between [a,b)
   */
  class RandomNumberDouble
  {
  public:
    RandomNumberDouble() = delete;

    RandomNumberDouble(const double _a, const double _b);

    double
    generate();

  private:
    double a, b;

    std::uniform_real_distribution<double> uniform_distribution;

    std::uint64_t timeSeed;

    std::seed_seq seed_sequence;

    std::mt19937_64 rng;
  };


  RandomNumberDouble::RandomNumberDouble(const double _a, const double _b)
    : a(_a)
    , b(_b)
    , uniform_distribution(a, b)
    , timeSeed(
        std::chrono::high_resolution_clock::now().time_since_epoch().count())
    , seed_sequence(
        {std::uint32_t(timeSeed & 0xffffffff), std::uint32_t(timeSeed >> 32)})
  {
    rng.seed(seed_sequence);
  }

  double
  RandomNumberDouble::generate()
  {
    return uniform_distribution(rng);
  }


  void
  fill_vector_randomly(Vector<double> &this_vector, double a = 0, double b = 1)
  {
    RandomNumberDouble random_double_generator(a, b);

    for (auto &v : this_vector)
      {
        v = random_double_generator.generate();
      }
  }
} // namespace FEConforimityTest
