// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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


// test deviator

#include <deal.II/base/symmetric_tensor.h>

#include "../tests.h"


template <int dim>
void
test()
{
  deallog << "dim=" << dim << std::endl;

  SymmetricTensor<2, dim> t;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      t[i][j] = (1. + (i + 1) * (j * 2));

  SymmetricTensor<2, dim> x = deviator(t);

  deallog << "x=" << std::endl;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      deallog << i << ' ' << j << ' ' << x[i][j] << std::endl;

  // the difference between t and x is a
  // diagonal tensor with constant elements
  // proportional to the trace of t
  t -= x;
  deallog << "t-x=" << std::endl;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      deallog << i << ' ' << j << ' ' << t[i][j] << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
