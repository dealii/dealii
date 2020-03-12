// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2018 by the deal.II authors
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


// check contract(Tensor<1,dim>,Tensor<2,dim>,Tensor<1,dim>)

#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

template <int dim>
void
test_select(double f1, double f2)
{
  Tensor<2, dim> t;
  unsigned int   k = 0;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      t[i][j] = ++k;

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      {
        Tensor<1, dim> u, v;
        u[i] = f1;
        v[j] = f2;
        deallog << '\t' << contract3(u, t, v);
      }
  deallog << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test_select<2>(1., 1.);
  test_select<3>(1., 1.);
  test_select<4>(1., 1.);

  test_select<2>(2., 3.);
  test_select<3>(2., 3.);
  test_select<4>(2., 3.);
}
