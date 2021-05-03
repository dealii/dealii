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


  // check product with Tensor<1,dim>
  check<Tensor<1, 2, Sdouble>, Tensor<1, 2, Sdouble>, Sdouble>();
  check<Tensor<1, 2, SSdouble>, Tensor<1, 2, SSdouble>, SSdouble>();

  Tensor<1, 2, SSdouble> t1;
  Tensor<1, 2, SSdouble> t2;
  SSdouble               a(2, 0, 7.0);
  SSdouble               b(2, 1, 3.0);
  SSdouble               c;
  a.val() = Sdouble(2, 0, 7.0);
  b.val() = Sdouble(2, 1, 3.0);

  for (unsigned int i = 0; i < 2; ++i)
    {
      t1[i] = 2. * a + i;
      t2[i] = 3. * b - i;
    }
  const Tensor<1, 2, SSdouble> t3 = t2;
  t1 *                         t2;

  t2 += a * t1;

  c = 0;
  c += (a * b + 0.5 * (t3 * t3 + t1 * t2)) * 0.3;



  deallog << "OK" << std::endl;
}
