// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/lac/vector.h>

#include "../tests.h"



const unsigned int N           = 10;
unsigned int       check_point = 0;



template <typename number>
void
print(const Vector<number> &v)
{
  //  deallog << "Check point " << check_point << std::endl;
  //  check_point++;

  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << v(i) << '\t';
  deallog << std::endl;
}



template <typename number1, typename number2>
void
check_vectors(Vector<number1> &d1, Vector<number2> &d2)
{
  deallog << "Fill & Swap" << std::endl;
  Vector<number1> d3(d1.size());
  print(d3);

  for (unsigned int i = 0; i < N; ++i)
    {
      d1(i) = 2. * i;
      d2(i) = .5 * d1(i) * d1(i);
      d3(i) = 2. - .5 * i;
    };

  print(d1);
  print(d2);
  print(d3);

  swap(d1, d2);
  print(d1);

  d1 = d2;
  print(d1);

  d1 = 2.5;
  print(d1);

  // initialize with iterators
  number1         array[] = {0.0, 1.1, 2.2, 3.3};
  Vector<number1> d4(&array[0], &array[4]);
  print(d4);

  deallog << "Extract number" << std::endl;
  // Each line should contain two equal numbers
  double sum = 0.;
  for (unsigned int i = 0; i < N; ++i)
    sum += 4. * i - i * i;
  deallog << d3 * d2 << '\t' << sum << std::endl;

  sum = 0.;
  for (unsigned int i = 0; i < N; ++i)
    sum += 4. * i * i;
  deallog << d2.norm_sqr() << '\t' << sum << std::endl;

  sum = std::sqrt(sum);
  deallog << d2.l2_norm() << '\t' << sum << std::endl;

  sum = 0.;
  for (unsigned int i = 0; i < N; ++i)
    sum += (2. - .5 * i) / N;
  deallog << d3.mean_value() << '\t' << sum << std::endl;

  sum = 0.;
  for (unsigned int i = 0; i < N; ++i)
    sum += std::fabs(2. - .5 * i);
  deallog << d3.l1_norm() << '\t' << sum << std::endl;

  sum = 0.;
  for (unsigned int i = 0; i < N; ++i)
    {
      double t = std::fabs(2. - .5 * i);
      if (t > sum)
        sum = t;
    }
  deallog << d3.linfty_norm() << '\t' << sum << std::endl;

  deallog << "add & sub" << std::endl;

  d1 += d2;
  print(d1);

  d2 -= d1;
  print(d2);

  d1.add(1.5);
  print(d1);

  d1.add(2, d3);
  print(d1);

  d1.add(2., d2, .5, d3);
  print(d1);

  deallog << "sadd & scale" << std::endl;

  d2.sadd(2., d1);
  print(d2);

  d2.sadd(2., .5, d1);
  print(d2);

  d1 *= 4.;
  print(d1);

  deallog << "equ" << std::endl;

  d2.equ(.25, d1);
  print(d2);
}


int
main()
{
  initlog();
  deallog << std::fixed;
  deallog << std::setprecision(2);

  Vector<double>      d1(N), d2(N);
  Vector<float>       f1(N), f2(N);
  Vector<long double> l1(N), l2(N);

  // cross-tests with double/float
  // vectors at the same time are not
  // supported at present,
  // unfortunately, as many functions
  // don't accept other data types as
  // arguments
  check_vectors(d1, d2);
  check_vectors(f1, f2);
  check_vectors(l1, l2);
}
