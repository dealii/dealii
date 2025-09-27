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


// test the determinant code

#include <deal.II/base/tensor.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Tensor<2, dim> t;

  // choose the same symmetric tensor
  // as in symmetric_tensor_09
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      t[i][j] = t[j][i] = (1. + (i + 1) * (j * 2));

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      deallog << i << ' ' << j << ' ' << t[i][j] << std::endl;

  deallog << determinant(t) << std::endl;
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
