// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// test for arithmetic operations on VectorizedArray

#include "../tests.h"
#include <iomanip>
#include <limits>

#include <deal.II/base/vectorization.h>


template <typename Number>
void test ()
{
  // since the number of array elements is system dependent, it is not a good
  // idea to print them to an output file. Instead, check the values manually
  VectorizedArray<Number> a, b, c;
  const unsigned int n_vectors = VectorizedArray<Number>::n_array_elements;
  a = Number(2.);
  b = Number(-1.);
  for (unsigned int i=0; i<n_vectors; ++i)
    c[i] = Number(i);

  deallog << "Addition: ";
  VectorizedArray<Number> d = a + b;
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (d[i] == 1, ExcInternalError());
  deallog << "OK" << std::endl
          << "Subtraction: ";
  VectorizedArray<Number> e = d - b;
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (e[i] == a[i], ExcInternalError());
  deallog << "OK" << std::endl
          << "Multiplication: ";
  d = a * c;
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (d[i] == a[i] * c[i], ExcInternalError());
  deallog << "OK" << std::endl
          << "Division: ";
  e = d / a;
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (e[i] == c[i], ExcInternalError());
  deallog << "OK" << std::endl
          << "Multiplication scalar: ";
  a = Number(2.) * a;
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (a[i] == 4., ExcInternalError());
  deallog << "OK" << std::endl
          << "Division scalar left: ";
  d = Number(1.) / a;
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (d[i] == 0.25, ExcInternalError());
  deallog << "OK" << std::endl
          << "Division scalar right: ";
  e = d / Number(0.25);
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (e[i] == 1, ExcInternalError());
  deallog << "OK" << std::endl
          << "Unary operator -: ";
  d = -c;
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (d[i] == -(Number)i, ExcInternalError());
  deallog << "OK" << std::endl
          << "Unary operator +: ";
  d = c;
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (d[i] == i, ExcInternalError());


  deallog << "OK" << std::endl
          << "Square root: ";
  d = std::sqrt(c);
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (std::fabs(d[i]-std::sqrt(Number(i)))<
            std::numeric_limits<Number>::epsilon(),
            ExcInternalError());

  deallog << "OK" << std::endl
          << "Absolute value: ";
  d = -c;
  d = std::abs(d);
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (d[i] == Number(i), ExcInternalError());
  deallog << "OK" << std::endl
          << "Maximum value: ";
  d = std::max(a, c);
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (d[i] == std::max(a[i], c[i]), ExcInternalError());
  deallog << "OK" << std::endl
          << "Minimum value: ";
  d = std::min(Number(0.5) * a + b, c);
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (d[i] == std::min(Number(0.5 * a[i] + b[i]), c[i]),
                 ExcInternalError());

  deallog << "OK" << std::endl
          << "Sine: ";
  e = std::sin(d);
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (std::fabs(e[i]-std::sin(d[i])) <
                 10.*std::numeric_limits<Number>::epsilon(),
                 ExcInternalError());
  deallog << "OK" << std::endl
          << "Cosine: ";
  e = std::cos(c);
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (std::fabs(e[i]-std::cos(c[i])) <
                 10.*std::numeric_limits<Number>::epsilon(),
                 ExcInternalError());
  deallog << "OK" << std::endl
          << "Tangent: ";
  d = std::tan(e);
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (std::fabs(d[i]-std::tan(e[i])) <
                 10.*std::numeric_limits<Number>::epsilon(),
                 ExcInternalError());
  deallog << "OK" << std::endl
          << "Exponential: ";
  d = std::exp(c-a);
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (std::fabs(d[i]-std::exp(c[i]-a[i])) <
                 10.*std::numeric_limits<Number>::epsilon(),
                 ExcInternalError());
  deallog << "OK" << std::endl
          << "Logarithm: ";
  e = std::log(d);
  for (unsigned int i=0; i<n_vectors; ++i)
    AssertThrow (std::fabs(e[i]-(c[i]-a[i])) <
                 10.*std::numeric_limits<Number>::epsilon(),
                 ExcInternalError());
  deallog << "OK" << std::endl;
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("double");
  test<double> ();
  deallog.pop();
  deallog.push("float");
  test<float> ();
  deallog.pop();

  // test long double: in that case, the default
  // path of VectorizedArray is taken no matter
  // what was done for double or float
  deallog.push("long double");
  test<long double> ();
  deallog.pop();
}
