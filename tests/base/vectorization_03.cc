// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2020 by the deal.II authors
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


// expression from step-48 that was computed in a wrong way for gcc-4.6.3 with
// vectorization


#include <deal.II/base/vectorization.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"


struct Evaluation
{
  VectorizedArray<double>
  get_value(const unsigned int index) const
  {
    return values[index];
  }

  void
  submit_value(const VectorizedArray<double> val, const unsigned int index)
  {
    if (is_cartesian)
      values[index] = val * cartesian_weight * jac_weight[index];
    else
      values[index] = val * general_weight[index];
  }

  bool                    is_cartesian;
  VectorizedArray<double> cartesian_weight;
  VectorizedArray<double> jac_weight[1];
  VectorizedArray<double> general_weight[1];
  VectorizedArray<double> values[1];
};


void
initialize(Evaluation &eval)
{
  eval.is_cartesian     = true;
  eval.cartesian_weight = random_value<double>();
  for (unsigned int i = 0; i < 4; ++i)
    eval.cartesian_weight =
      std::max(eval.cartesian_weight,
               eval.cartesian_weight * eval.cartesian_weight);
  eval.general_weight[0] = 0.2313342 * eval.cartesian_weight;
  eval.jac_weight[0]     = random_value<double>();
}


void
test()
{
  Evaluation current, old;
  initialize(current);
  initialize(old);
  VectorizedArray<double> weight;
  weight = random_value<double>();

  VectorizedArray<double> vec;
  for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v)
    vec[v] = random_value<double>();

  current.values[0] = vec;
  old.values[0] =
    vec * 1.112 - std::max(2. * vec - 1., VectorizedArray<double>());

  Vector<double> vector(200);
  vector = 1.2;

  VectorizedArray<double> cur = current.get_value(0);
  VectorizedArray<double> ol  = old.get_value(0);
  current.submit_value(2. * cur - ol - weight * std::sin(cur), 0);

  vector *= 2. * current.get_value(0)[0];

  double error = 0;
  for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v)
    error += std::abs(current.get_value(0)[v] / (current.cartesian_weight[v] *
                                                 current.jac_weight[0][0]) -
                      (2. * vec[v] - ol[v] - weight[v] * std::sin(vec[v])));
  deallog << "error: " << error << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(4);

  test();
}
