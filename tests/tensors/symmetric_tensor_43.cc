// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Verify that symmetric tensor operators+= and friends work with different
// number types

#include <deal.II/base/symmetric_tensor.h>

#include <complex>

#include "../tests.h"

template <int rank, int dim, typename NumberType>
void
fill_tensor(SymmetricTensor<rank, dim, NumberType> &t)
{
  for (unsigned int i = 0; i != t.n_independent_components; ++i)
    {
      t.access_raw_entry(i) = i + 1;
    }
}

template <int dim, int rank, typename NumberType, typename OtherNumberType>
void
test_operators()
{
  const OtherNumberType s1 = 2.0;
  const OtherNumberType s2 = 0.5;

  deallog << "operator *=" << std::endl;
  SymmetricTensor<rank, dim, NumberType> t1;
  fill_tensor(t1);
  t1 *= s1;
  deallog << t1 << std::endl;

  deallog << "operator /=" << std::endl;
  SymmetricTensor<rank, dim, NumberType> t2;
  fill_tensor(t2);
  t2 /= s2;
  deallog << t2 << std::endl;

  SymmetricTensor<rank, dim, OtherNumberType> t3;
  fill_tensor(t3);

  deallog << "operator +=" << std::endl;
  t1 += t3;
  deallog << t1 << std::endl;

  deallog << "operator -=" << std::endl;
  t2 -= t3;
  deallog << t2 << std::endl;
}

template <int rank, typename NumberType, typename OtherNumberType>
void
test()
{
  deallog << "dim=" << 1 << std::endl;
  test_operators<1, rank, NumberType, OtherNumberType>();
  deallog << "dim=" << 2 << std::endl;
  test_operators<2, rank, NumberType, OtherNumberType>();
  deallog << "dim=" << 3 << std::endl;
  test_operators<3, rank, NumberType, OtherNumberType>();
}


int
main()
{
  initlog();

  deallog << "rank 2, double" << std::endl;
  test<2, double, float>();
  deallog << "rank 2, std::complex<double>" << std::endl;
  test<2, std::complex<double>, std::complex<float>>();

  deallog << "rank 4, double" << std::endl;
  test<4, double, float>();
  deallog << "rank 4, std::complex<double>" << std::endl;
  test<4, std::complex<double>, std::complex<float>>();

  deallog << "OK" << std::endl;
}
