// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// test that the product between a tensor of rank 1 and a scalar
// results a result type as expected

#include <deal.II/base/template_constraints.h>
#include <deal.II/base/tensor.h>

#include <complex>
#include <typeinfo>

#include "../tests.h"


template <typename T, typename U, typename CompareType>
void
check()
{
  AssertThrow(typeid(T() * U()) == typeid(CompareType), ExcInternalError());
  AssertThrow(typeid(T() * U()) == typeid(CompareType), ExcInternalError());
}


int
main()
{
  initlog();

  // check product of scalars
  check<Tensor<1, 1, double>, double, Tensor<1, 1, double>>();
  check<double, Tensor<1, 1, double>, Tensor<1, 1, double>>();

  check<Tensor<1, 1, double>, float, Tensor<1, 1, double>>();
  check<float, Tensor<1, 1, double>, Tensor<1, 1, double>>();

  check<Tensor<1, 1, double>,
        std::complex<double>,
        Tensor<1, 1, std::complex<double>>>();
  check<std::complex<double>,
        Tensor<1, 1, double>,
        Tensor<1, 1, std::complex<double>>>();

  check<Tensor<1, 1, double>,
        std::complex<float>,
        Tensor<1, 1, std::complex<double>>>();
  check<std::complex<float>,
        Tensor<1, 1, double>,
        Tensor<1, 1, std::complex<double>>>();

  check<Tensor<1, 1, float>,
        std::complex<double>,
        Tensor<1, 1, std::complex<double>>>();
  check<std::complex<double>,
        Tensor<1, 1, float>,
        Tensor<1, 1, std::complex<double>>>();

  check<Tensor<1, 1, float>,
        std::complex<float>,
        Tensor<1, 1, std::complex<float>>>();
  check<std::complex<float>,
        Tensor<1, 1, float>,
        Tensor<1, 1, std::complex<float>>>();

  deallog << "OK" << std::endl;
}
