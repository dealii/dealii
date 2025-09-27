// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Tests LAPACKFullMatrix::transpose(...) const

#include <deal.II/lac/lapack_full_matrix.h>

#include <iostream>
#include <tuple>

#include "../tests.h"

#include "create_matrix.h"

template <typename T>
std::string
to_string(const T &t)
{
  std::ostringstream s;
  s << t;
  return s.str();
}

template <typename NumberType>
void
test(const unsigned int n = 3, const unsigned int k = 6)
{
  LAPACKFullMatrix<NumberType> A(n, k);
  LAPACKFullMatrix<NumberType> At(k, n);

  create_random(A);
  A.transpose(At);
  for (unsigned int i = 0; i < A.m(); ++i)
    for (unsigned int j = 0; j < A.n(); ++j)
      {
        const auto at = At(j, i);
        const auto a  = numbers::NumberTraits<NumberType>::conjugate(A(i, j));
        AssertThrow(at == a,
                    ExcMessage(to_string(a) + "!=" + to_string(at) + " for (" +
                               std::to_string(i) + "," + std::to_string(j) +
                               ")"));
      }

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  test<double>(11, 27);
  test<double>(15, 4);
  test<double>(10, 10);

  test<std::complex<double>>(11, 27);
  test<std::complex<double>>(15, 4);
  test<std::complex<double>>(10, 10);
}
