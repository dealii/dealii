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


// like _01 but for SymmetricTensor

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/template_constraints.h>

#include <complex>
#include <typeinfo>

#include "../tests.h"


template <typename T, typename U, typename CompareType>
void
check()
{
  AssertThrow(typeid(typename ProductType<T, U>::type) == typeid(CompareType),
              ExcInternalError());
  AssertThrow(typeid(typename ProductType<T, U>::type) == typeid(T() * U()),
              ExcInternalError());
}


int
main()
{
  initlog();

  // check product with SymmetricTensor<2,dim>
  check<SymmetricTensor<2, 2, double>, double, SymmetricTensor<2, 2, double>>();
  check<SymmetricTensor<2, 2, float>, double, SymmetricTensor<2, 2, double>>();
  check<double, SymmetricTensor<2, 2, float>, SymmetricTensor<2, 2, double>>();

  // check product with std::complex. rather annoyingly, there is no
  // product between std::complex<double> and float, or the other way
  // around, so stay within the same type system
  check<std::complex<double>, double, std::complex<double>>();
  check<std::complex<float>, float, std::complex<float>>();
  check<SymmetricTensor<2, 2>,
        std::complex<double>,
        SymmetricTensor<2, 2, std::complex<double>>>();
  check<std::complex<double>,
        SymmetricTensor<2, 2>,
        SymmetricTensor<2, 2, std::complex<double>>>();

  deallog << "OK" << std::endl;
}
