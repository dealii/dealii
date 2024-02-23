// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
