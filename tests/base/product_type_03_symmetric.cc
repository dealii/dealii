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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


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
