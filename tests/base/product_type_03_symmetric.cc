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


// like _03 but for SymmetricTensor

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/template_constraints.h>

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
  check<SymmetricTensor<2, 1, double>, double, SymmetricTensor<2, 1, double>>();
  check<double, SymmetricTensor<2, 1, double>, SymmetricTensor<2, 1, double>>();

  check<SymmetricTensor<2, 1, double>, float, SymmetricTensor<2, 1, double>>();
  check<float, SymmetricTensor<2, 1, double>, SymmetricTensor<2, 1, double>>();

  check<SymmetricTensor<2, 1, double>,
        std::complex<double>,
        SymmetricTensor<2, 1, std::complex<double>>>();
  check<std::complex<double>,
        SymmetricTensor<2, 1, double>,
        SymmetricTensor<2, 1, std::complex<double>>>();

  check<SymmetricTensor<2, 1, double>,
        std::complex<float>,
        SymmetricTensor<2, 1, std::complex<double>>>();
  check<std::complex<float>,
        SymmetricTensor<2, 1, double>,
        SymmetricTensor<2, 1, std::complex<double>>>();

  check<SymmetricTensor<2, 1, float>,
        std::complex<double>,
        SymmetricTensor<2, 1, std::complex<double>>>();
  check<std::complex<double>,
        SymmetricTensor<2, 1, float>,
        SymmetricTensor<2, 1, std::complex<double>>>();

  check<SymmetricTensor<2, 1, float>,
        std::complex<float>,
        SymmetricTensor<2, 1, std::complex<float>>>();
  check<std::complex<float>,
        SymmetricTensor<2, 1, float>,
        SymmetricTensor<2, 1, std::complex<float>>>();

  deallog << "OK" << std::endl;
}
