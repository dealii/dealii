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


// Verify that symmetric tensor construction and initialization
// works as expected for non-trivial number types

#include <deal.II/base/symmetric_tensor.h>

#include <complex>

#include "../tests.h"

template <int dim, int rank, typename NumberType>
void
test_constructor()
{
  const SymmetricTensor<rank, dim, NumberType> t1;
}

template <int dim, int rank, typename NumberType, typename OtherNumberType>
void
test_initializer()
{
  const SymmetricTensor<rank, dim, OtherNumberType> t1;
  const SymmetricTensor<rank, dim, NumberType>      t2(t1);
}

template <int rank, typename NumberType, typename OtherNumberType>
void
test()
{
  deallog << "constructor" << std::endl;
  test_constructor<1, rank, NumberType>();
  test_constructor<2, rank, NumberType>();
  test_constructor<3, rank, NumberType>();

  deallog << "initializer" << std::endl;
  test_initializer<1, rank, NumberType, OtherNumberType>();
  test_initializer<2, rank, NumberType, OtherNumberType>();
  test_initializer<3, rank, NumberType, OtherNumberType>();
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
