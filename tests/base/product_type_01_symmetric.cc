// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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
