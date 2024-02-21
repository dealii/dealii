// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// checks that the l1, l2, lp norm is computed correctly for deal.II vectors
// (check the summation algorithm), including an accuracy test (should not
// lose more than 1 decimal also for 200000 vector entries)

#include <deal.II/base/numbers.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



template <typename number>
void
check_norms()
{
  const number acc  = 1e1 * std::numeric_limits<number>::epsilon();
  unsigned int skip = 73;
  for (unsigned int size = 1; size < 200000; size += skip)
    {
      // test correctness
      if (size > 10000)
        skip += 17;
      Vector<number> vec(size);
      for (unsigned int i = 0; i < size; ++i)
        vec(i) = i + 1;

      const number l1_norm = vec.l1_norm();
      AssertThrow(std::abs(l1_norm - 0.5 * size * (size + 1)) <
                    acc * 0.5 * size * (size + 1),
                  ExcInternalError());

      // test accuracy of summation
      constexpr long double value = numbers::PI;
      vec                         = (number)value;
      const number l1_norma       = vec.l1_norm();
      AssertThrow(std::abs(l1_norma - value * size) < acc * size * value,
                  ExcInternalError());
      const number l2_norma = vec.l2_norm();
      AssertThrow(std::abs(l2_norma - value * std::sqrt((number)size)) <
                    acc * std::sqrt(size) * value,
                  ExcInternalError());
      const number lp_norma = vec.lp_norm(3.);
      AssertThrow(
        std::abs(lp_norma - value * std::pow(static_cast<number>(size),
                                             static_cast<number>(1. / 3.))) <
          std::max(acc, number(1e-15)) *
            std::pow(static_cast<number>(size), static_cast<number>(1. / 3.)) *
            value,
        ExcInternalError());
    }
}


template <typename number>
void
check_complex_norms()
{
  const number acc  = 1e2 * std::numeric_limits<number>::epsilon();
  unsigned int skip = 73;
  for (unsigned int size = 1; size < 100000; size += skip)
    {
      // test correctness
      if (size > 10000)
        skip += 17;
      Vector<std::complex<number>> vec(size);
      long double                  sum = 0.;
      for (unsigned int i = 0; i < size; ++i)
        {
          vec(i) = std::complex<number>(i + 1, i + 2);
          sum += std::sqrt((long double)(i + 1) * (i + 1) +
                           (long double)(i + 2) * (i + 2));
        }

      const number l1_norm = vec.l1_norm();
      AssertThrow(std::abs(l1_norm - sum) < acc * sum, ExcInternalError());

      // test accuracy of summation
      constexpr std::complex<long double> value(numbers::PI, 0.1);
      vec                   = std::complex<number>(value);
      const number l1_norma = vec.l1_norm();
      AssertThrow(std::abs(l1_norma - std::abs(value) * size) <
                    acc * size * std::abs(value),
                  ExcInternalError());
      const number l2_norma = vec.l2_norm();
      AssertThrow(std::abs(l2_norma -
                           std::abs(value) * std::sqrt((number)size)) <
                    acc * std::sqrt((number)size) * std::abs(value),
                  ExcInternalError());
      const number lp_norma = vec.lp_norm(3.);
      AssertThrow(
        std::abs(lp_norma -
                 std::abs(value) * std::pow(static_cast<number>(size),
                                            static_cast<number>(1. / 3.))) <
          acc *
            std::pow(static_cast<number>(size), static_cast<number>(1. / 3.)) *
            std::abs(value),
        ExcInternalError());
    }
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);

  check_norms<float>();
  check_norms<double>();
  check_norms<long double>();
#ifdef DEAL_II_WITH_COMPLEX_VALUES
  check_complex_norms<double>();
  check_complex_norms<float>();
#endif
  deallog << "OK" << std::endl;
}
