// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// a stress test using repeated multiplication of tensors

#include <deal.II/base/tensor.h>

#include "../tests.h"



template <int dim>
void
test()
{
  const double   lambda = 1.5, mu = 1.7;
  Tensor<2, dim> t;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      t[i][j] = (1. + (i + lambda) * (mu + 13));

  Tensor<2, dim> a;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      a[i][j] = (1. + (i + 1) * (j + 1));

  // stress test the whole thing many
  // times. normalize in each step to
  // make sure the result remains
  // representable in floating point
  // arithmetic. essentially, this
  // invokes the power method to
  // compute the largest eigenvector
  // (eigentensor in this case)
  for (unsigned int i = 0; i < 1000000; ++i)
    {
      a = t * a;
      a /= a.norm();
    }

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      deallog << i << ' ' << j << ' ' << a[i][j] << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
