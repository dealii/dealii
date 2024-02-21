// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test the SymbolicFunction class

#include <deal.II/base/symbolic_function.h>

#include <map>

#include "../tests.h"

using namespace Differentiation::SD;

template <int dim>
Functions::SymbolicFunction<dim>
test(const std::string &expression)
{
  auto x  = Functions::SymbolicFunction<dim>::get_default_coordinate_symbols();
  auto t  = make_symbol("t");
  auto f  = Expression(expression, true);
  auto df = differentiate(f, x);
  auto Hf = differentiate(df, x);

  // Create a Function with a single component and an evaluation point
  Functions::SymbolicFunction<dim> fun({f}, x);
  Point<dim>                       p;
  for (unsigned int i = 0; i < dim; ++i)
    p[i] = i + 1.0;

  // Output all symbolic stuff
  deallog << "=========================================================="
          << std::endl
          << "dim = " << dim << ", Symengine" << std::endl
          << "x: " << x << ", f: " << f << ", df: " << df << ", H:" << Hf
          << std::endl;

  // Output the function and its evaluation
  deallog << "SymbolicFunction<dim>: " << fun << std::endl
          << "p: " << p << ", f(p): " << fun.value(p)
          << ", grad(f)(p): " << fun.gradient(p)
          << ", laplacian(f)(p): " << fun.laplacian(p)
          << ", H(f)(p): " << fun.hessian(p) << std::endl;

  // Output its time derivative
  auto fun_t = fun.time_derivative();
  deallog << "Time derivative of SymbolicFunction<dim>: " << fun_t << std::endl
          << "p: " << p << ", f_t(p): " << fun_t.value(p)
          << ", grad(f_t)(p): " << fun_t.gradient(p)
          << ", laplacian(f_t)(p): " << fun_t.laplacian(p)
          << ", H(f_t)(p): " << fun_t.hessian(p) << std::endl;

  return fun;
}

int
main()
{
  initlog();

  test<1>("x**2 + t*x");
  test<2>("x+y+cos(exp(y*x*t))");
  test<3>("tan(z*t**2) + atan(y/x)");
}
