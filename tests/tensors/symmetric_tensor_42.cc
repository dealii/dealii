// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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
