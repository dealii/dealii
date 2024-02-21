// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test that the product between a tensor of rank 2 and a scalar
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
  check<Tensor<2, 1, double>, double, Tensor<2, 1, double>>();
  check<double, Tensor<2, 1, double>, Tensor<2, 1, double>>();

  check<Tensor<2, 1, double>, float, Tensor<2, 1, double>>();
  check<float, Tensor<2, 1, double>, Tensor<2, 1, double>>();

  check<Tensor<2, 1, double>,
        std::complex<double>,
        Tensor<2, 1, std::complex<double>>>();
  check<std::complex<double>,
        Tensor<2, 1, double>,
        Tensor<2, 1, std::complex<double>>>();

  check<Tensor<2, 1, double>,
        std::complex<float>,
        Tensor<2, 1, std::complex<double>>>();
  check<std::complex<float>,
        Tensor<2, 1, double>,
        Tensor<2, 1, std::complex<double>>>();

  check<Tensor<2, 1, float>,
        std::complex<double>,
        Tensor<2, 1, std::complex<double>>>();
  check<std::complex<double>,
        Tensor<2, 1, float>,
        Tensor<2, 1, std::complex<double>>>();

  check<Tensor<2, 1, float>,
        std::complex<float>,
        Tensor<2, 1, std::complex<float>>>();
  check<std::complex<float>,
        Tensor<2, 1, float>,
        Tensor<2, 1, std::complex<float>>>();

  deallog << "OK" << std::endl;
}
