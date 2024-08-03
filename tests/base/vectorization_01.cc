// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for arithmetic operations on VectorizedArray

#include <deal.II/base/vectorization.h>

#include <limits>

#include "../tests.h"


template <typename Number>
void
test()
{
  // since the number of array elements is system dependent, it is not a good
  // idea to print them to an output file. Instead, check the values manually
  VectorizedArray<Number> a, b, c;
  const unsigned int      n_vectors = VectorizedArray<Number>::size();
  a                                 = Number(2.);
  b                                 = Number(-1.);
  for (unsigned int i = 0; i < n_vectors; ++i)
    c[i] = Number(i);

  AssertDimension(n_vectors, a.size());
  AssertDimension(a.size(), sizeof(a) / sizeof(Number));
  AssertDimension(VectorizedArray<Number>::size(), sizeof(a) / sizeof(Number));

  deallog << "Addition: ";
  VectorizedArray<Number> d = a + b;
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(d[i] == 1, ExcInternalError());
  deallog << "OK" << std::endl << "Subtraction: ";
  VectorizedArray<Number> e = d - b;
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(e[i] == a[i], ExcInternalError());
  deallog << "OK" << std::endl << "Multiplication: ";
  d = a * c;
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(d[i] == a[i] * c[i], ExcInternalError());
  deallog << "OK" << std::endl << "Division: ";
  e = d / a;
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(e[i] == c[i], ExcInternalError());
  deallog << "OK" << std::endl << "Multiplication scalar: ";
  a = Number(2.) * a;
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(a[i] == 4., ExcInternalError());
  deallog << "OK" << std::endl << "Division scalar left: ";
  d = Number(1.) / a;
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(d[i] == 0.25, ExcInternalError());
  deallog << "OK" << std::endl << "Division scalar right: ";
  e = d / Number(0.25);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(e[i] == 1, ExcInternalError());
  deallog << "OK" << std::endl << "Unary operator -: ";
  d = -c;
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(d[i] == -(Number)i, ExcInternalError());
  deallog << "OK" << std::endl << "Unary operator +: ";
  d = c;
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(d[i] == i, ExcInternalError());


  deallog << "OK" << std::endl << "Square root: ";
  d = std::sqrt(c);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(std::fabs(d[i] - std::sqrt(Number(i))) <
                  std::numeric_limits<Number>::epsilon(),
                ExcInternalError());

  deallog << "OK" << std::endl << "Absolute value: ";
  d = -c;
  d = std::abs(d);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(d[i] == Number(i), ExcInternalError());
  deallog << "OK" << std::endl << "Maximum value: ";
  d = std::max(a, c);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(d[i] == std::max(a[i], c[i]), ExcInternalError());
  deallog << "OK" << std::endl << "Minimum value: ";
  d = std::min(Number(0.5) * a + b, c);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(d[i] == std::min(Number(0.5 * a[i] + b[i]), c[i]),
                ExcInternalError());

  deallog << "OK" << std::endl << "Sine: ";
  e = std::sin(d);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(std::fabs(e[i] - std::sin(d[i])) <
                  10. * std::numeric_limits<Number>::epsilon(),
                ExcInternalError());
  deallog << "OK" << std::endl << "Arc sine: ";
  d = std::asin(e);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(std::fabs(d[i] - std::asin(e[i])) <
                  10. * std::numeric_limits<Number>::epsilon(),
                ExcInternalError());
  deallog << "OK" << std::endl << "Cosine: ";
  e = std::cos(c);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(std::fabs(e[i] - std::cos(c[i])) <
                  10. * std::numeric_limits<Number>::epsilon(),
                ExcInternalError());
  deallog << "OK" << std::endl << "Arc cosine: ";
  c = std::acos(e);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(std::fabs(c[i] - std::acos(e[i])) <
                  10. * std::numeric_limits<Number>::epsilon(),
                ExcInternalError());
  deallog << "OK" << std::endl << "Tangent: ";
  d = std::tan(e);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(std::fabs(d[i] - std::tan(e[i])) <
                  10. * std::numeric_limits<Number>::epsilon(),
                ExcInternalError());
  deallog << "OK" << std::endl << "Arc tangent: ";
  d = std::atan(e);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(std::fabs(d[i] - std::atan(e[i])) <
                  10. * std::numeric_limits<Number>::epsilon(),
                ExcInternalError());
  deallog << "OK" << std::endl << "Hyperbolic cosine: ";
  e = std::cosh(c);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(std::fabs(e[i] - std::cosh(c[i])) <
                  10. * std::numeric_limits<Number>::epsilon(),
                ExcInternalError());
  deallog << "OK" << std::endl << "Area hyperbolic cosine: ";
  c = std::acosh(e);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(std::fabs(c[i] - std::acosh(e[i])) <
                  10. * std::numeric_limits<Number>::epsilon(),
                ExcInternalError());
  deallog << "OK" << std::endl << "Hyperbolic sine: ";
  e = std::sinh(d);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(std::fabs(e[i] - std::sinh(d[i])) <
                  10. * std::numeric_limits<Number>::epsilon(),
                ExcInternalError());
  deallog << "OK" << std::endl << "Area hyperbolic sine: ";
  d = std::asinh(e);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(std::fabs(d[i] - std::asinh(e[i])) <
                  10. * std::numeric_limits<Number>::epsilon(),
                ExcInternalError());
  deallog << "OK" << std::endl << "Hyperbolic tangent: ";
  d = std::tanh(e);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(std::fabs(d[i] - std::tanh(e[i])) <
                  10. * std::numeric_limits<Number>::epsilon(),
                ExcInternalError());
  deallog << "OK" << std::endl << "Area hyperbolic tangent: ";
  e = std::atanh(d);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(std::fabs(e[i] - std::atanh(d[i])) <
                  10. * std::numeric_limits<Number>::epsilon(),
                ExcInternalError());
  deallog << "OK" << std::endl << "Exponential: ";
  d = std::exp(c - a);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(std::fabs(d[i] - std::exp(c[i] - a[i])) <
                  10. * std::numeric_limits<Number>::epsilon(),
                ExcInternalError());
  deallog << "OK" << std::endl << "Logarithm: ";
  e = std::log(d);
  for (unsigned int i = 0; i < n_vectors; ++i)
    AssertThrow(std::fabs(e[i] - (c[i] - a[i])) <
                  10. * std::numeric_limits<Number>::epsilon(),
                ExcInternalError());
  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  deallog.push("double");
  test<double>();
  deallog.pop();
  deallog.push("float");
  test<float>();
  deallog.pop();

  // test long double: in that case, the default
  // path of VectorizedArray is taken no matter
  // what was done for double or float
  deallog.push("long double");
  test<long double>();
  deallog.pop();
}
