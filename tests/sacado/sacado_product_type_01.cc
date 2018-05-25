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
  typedef Sacado::Fad::DFad<double> Sdouble;
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
