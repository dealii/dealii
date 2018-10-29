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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// test scalar_product between tensors and symmetric tensors with Sacado


#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/differentiation/ad/sacado_product_types.h>

#include "../tests.h"


int
main()
{
  typedef Sacado::Fad::DFad<double>  Sdouble;
  typedef Sacado::Fad::DFad<Sdouble> SSdouble;
  initlog();


  // check product with Tensor<2,dim>
  Tensor<2, 2, SSdouble>        t;
  SymmetricTensor<2, 2, double> st;
  SSdouble                      a(2, 0, 7.0);
  SSdouble                      b(2, 1, 3.0);
  a.val() = Sdouble(2, 0, 7.0);
  b.val() = Sdouble(2, 1, 3.0);

  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      {
        t[i][j]  = 2. * a + i * j * b;
        st[i][j] = (1. + (i + 1) * (j * 2));
      }

  deallog << scalar_product(t, st) << std::endl;
  deallog << scalar_product(st, t) << std::endl;
}
