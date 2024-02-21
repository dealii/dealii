// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check FunctionFromFunctionObjects implementation

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/vector.h>

#include <boost/core/demangle.hpp>

#include <string>
#include <typeinfo>
#include <vector>

#include "../tests.h"

// Test the lambda function object.

int
main()
{
  initlog();

  // Scalar function:
  auto val1 = [](const Point<2> &p) { return p.square(); };
  auto val2 = [](const Point<2> &p) { return 2 * p[0]; };

  auto grad1 = [](const Point<2> &p) { return 2 * Tensor<1, 2>(p); };
  auto grad2 = [](const Point<2> &) { return 2 * Tensor<1, 2>({1, 0}); };

  FunctionFromFunctionObjects<2> fun1({val1});
  FunctionFromFunctionObjects<2> fun2({val1, val2});

  FunctionFromFunctionObjects<2> fun3({val1}, {grad1});
  FunctionFromFunctionObjects<2> fun4({val1, val2}, {grad1, grad2});


  Point<2> p(1, 2);

  deallog << fun1.value(p) << std::endl;

  Vector<double> v(2);
  fun2.vector_value(p, v);
  deallog << v(0) << ' ' << v(1) << std::endl;

  deallog << fun3.gradient(p) << std::endl;

  std::vector<Tensor<1, 2>> g(2);
  fun4.vector_gradient(p, g);
  deallog << g[0] << ", " << g[1] << std::endl;
}
