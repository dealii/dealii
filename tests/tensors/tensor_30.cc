// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Additional tests for tensor contraction
//
// A mysteriously failing pull request suggested that the contraction
// of a rank-3 tensor and a rank-1 tensor might be using the wrong
// indices. Check this.

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Tensor<3, dim> a;
  Tensor<1, dim> b;

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        a[i][j][k] = i + 2 * j + 3 * k;

  for (unsigned int k = 0; k < dim; ++k)
    b[k] = k;

  // Compute
  //   c(i,j) = sum_k a(i,j,k)*b(k)
  //          = sum_k (i + 2j + 3k)k
  //          = (i + 2j) sum k  + 3 sum k^2
  // which in 2d equates to
  //          = (i + 2j) + 3
  // and in 3d equates to
  //          = 3(i + 2j) + 3*5
  //
  // In other words, in 2d this is the 2x2 matrix
  //    [[3,4],[4,6]]
  // and in 3d is the 3x3 matrix
  //    [[15,21,27],[18,24,30],[21,27,33]]
  Tensor<2, dim> c = a * b;

  deallog << c << std::endl;
}



int
main(int argc, char *argv[])
{
  initlog();

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;

  return 0;
}
