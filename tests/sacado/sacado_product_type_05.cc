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

#include <typeinfo>

#include <deal.II/base/tensor.h>
#include <deal.II/differentiation/ad/sacado_product_types.h>

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
  typedef Sacado::Fad::DFad<double>  Sdouble;
  typedef Sacado::Fad::DFad<Sdouble> SSdouble;
  initlog();

  // check product with Tensor<2,dim>
  check<Tensor<1, 2, Sdouble>, Tensor<1, 2, Sdouble>, Sdouble>();
  check<Tensor<1, 2, SSdouble>, Tensor<1, 2, SSdouble>, SSdouble>();
  check<Tensor<2, 2, Sdouble>, Tensor<2, 2, Sdouble>, Tensor<2, 2, Sdouble>>();
  check<Tensor<2, 2, SSdouble>,
        Tensor<2, 2, SSdouble>,
        Tensor<2, 2, SSdouble>>();

  Tensor<2, 2, SSdouble> t1;
  Tensor<2, 2, SSdouble> t2;
  SSdouble               a(2, 0, 7.0);
  SSdouble               b(2, 1, 3.0);
  SSdouble               c;
  a.val() = Sdouble(2, 0, 7.0);
  b.val() = Sdouble(2, 1, 3.0);

  for(unsigned int i = 0; i < 2; ++i)
    for(unsigned int j = 0; j < 2; ++j)
      {
        t1[i][j] = 2. * a + i * j;
        t2[i][j] = 3. * b - i + j;
      }
  t1* t2;

  t2 += a * t1;

  deallog << "OK" << std::endl;
}
