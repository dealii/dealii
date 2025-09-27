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
  using Sdouble  = Sacado::Fad::DFad<double>;
  using SSdouble = Sacado::Fad::DFad<Sdouble>;
  initlog();


  // check product of scalars
  check<Sdouble, SSdouble, SSdouble>();
  check<SSdouble, Sdouble, SSdouble>();

  // check product with Tensor<1,dim>
  check<Tensor<1, 2, Sdouble>, SSdouble, Tensor<1, 2, SSdouble>>();
  check<Tensor<1, 2, SSdouble>, Sdouble, Tensor<1, 2, SSdouble>>();

  check<SSdouble, Tensor<1, 2, Sdouble>, Tensor<1, 2, SSdouble>>();
  check<Sdouble, Tensor<1, 2, SSdouble>, Tensor<1, 2, SSdouble>>();

  // check product with Tensor<2,dim> (which is a different class)
  check<Tensor<2, 2, Sdouble>, SSdouble, Tensor<2, 2, SSdouble>>();
  check<Tensor<2, 2, SSdouble>, Sdouble, Tensor<2, 2, SSdouble>>();

  check<SSdouble, Tensor<2, 2, Sdouble>, Tensor<2, 2, SSdouble>>();
  check<Sdouble, Tensor<2, 2, SSdouble>, Tensor<2, 2, SSdouble>>();

  deallog << "OK" << std::endl;
}
