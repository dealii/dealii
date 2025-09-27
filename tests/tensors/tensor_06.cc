// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check Tensor<1,dim>::operator*(Tensor<1,dim>)

#include <deal.II/base/tensor.h>

#include "../tests.h"

template <int dim>
void
test_tensor()
{
  Tensor<1, dim> t1, t2;
  for (unsigned int i = 0; i < dim; ++i)
    {
      t1[i] = i + 1;
      t2[i] = 4. * i - 10.;
    }
  double res = t1 * t2;
  deallog << "dim = " << dim << ": " << res << std::endl;
}

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test_tensor<1>();
  test_tensor<2>();
  test_tensor<3>();
  test_tensor<4>();
  test_tensor<7>();
  deallog << "OK" << std::endl;
}
