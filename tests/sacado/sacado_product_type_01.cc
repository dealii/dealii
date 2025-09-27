// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test ProductType with sacado

#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/ad/sacado_product_types.h>

#include <typeinfo>

#include "../tests.h"

template <typename T, typename U, typename CompareType>
void
check()
{
  AssertThrow(typeid(typename ProductType<T, U>::type) == typeid(CompareType),
              ExcInternalError());
}


int
main()
{
  using Sdouble = Sacado::Fad::DFad<double>;
  initlog();


  // check product of scalars
  check<double, Sdouble, Sdouble>();
  check<Sdouble, double, Sdouble>();

  check<float, Sdouble, Sdouble>();
  check<Sdouble, float, Sdouble>();

  check<int, Sdouble, Sdouble>();
  check<Sdouble, int, Sdouble>();

  // check product with Tensor<1,dim>
  check<Tensor<1, 2, double>, Sdouble, Tensor<1, 2, Sdouble>>();
  check<Tensor<1, 2, Sdouble>, double, Tensor<1, 2, Sdouble>>();

  check<Sdouble, Tensor<1, 2, double>, Tensor<1, 2, Sdouble>>();
  check<double, Tensor<1, 2, Sdouble>, Tensor<1, 2, Sdouble>>();

  check<Tensor<1, 2, float>, Sdouble, Tensor<1, 2, Sdouble>>();
  check<Tensor<1, 2, Sdouble>, float, Tensor<1, 2, Sdouble>>();

  check<Sdouble, Tensor<1, 2, float>, Tensor<1, 2, Sdouble>>();
  check<float, Tensor<1, 2, Sdouble>, Tensor<1, 2, Sdouble>>();

  check<Tensor<1, 2, int>, Sdouble, Tensor<1, 2, Sdouble>>();
  check<Tensor<1, 2, Sdouble>, int, Tensor<1, 2, Sdouble>>();

  check<Sdouble, Tensor<1, 2, int>, Tensor<1, 2, Sdouble>>();
  check<int, Tensor<1, 2, Sdouble>, Tensor<1, 2, Sdouble>>();

  // check product with Tensor<2,dim> (which is a different class)
  check<Tensor<2, 2, double>, Sdouble, Tensor<2, 2, Sdouble>>();
  check<Tensor<2, 2, Sdouble>, double, Tensor<2, 2, Sdouble>>();

  check<Sdouble, Tensor<2, 2, double>, Tensor<2, 2, Sdouble>>();
  check<double, Tensor<2, 2, Sdouble>, Tensor<2, 2, Sdouble>>();

  check<Tensor<2, 2, float>, Sdouble, Tensor<2, 2, Sdouble>>();
  check<Tensor<2, 2, Sdouble>, float, Tensor<2, 2, Sdouble>>();

  check<Sdouble, Tensor<2, 2, float>, Tensor<2, 2, Sdouble>>();
  check<float, Tensor<2, 2, Sdouble>, Tensor<2, 2, Sdouble>>();

  check<Tensor<2, 2, int>, Sdouble, Tensor<2, 2, Sdouble>>();
  check<Tensor<2, 2, Sdouble>, int, Tensor<2, 2, Sdouble>>();

  check<Sdouble, Tensor<2, 2, int>, Tensor<2, 2, Sdouble>>();
  check<int, Tensor<2, 2, Sdouble>, Tensor<2, 2, Sdouble>>();

  deallog << "OK" << std::endl;
}
