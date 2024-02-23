// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test outer_product

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

  // test 1: check trace-like operator
  {
    const SymmetricTensor<4, dim> T =
      outer_product<dim>(unit_symmetric_tensor<dim>(),
                         unit_symmetric_tensor<dim>());

    // T*t should yield a diagonal tensor
    // where the diagonal elements are the
    // traces of t
    SymmetricTensor<2, dim> x = T * t;
    AssertThrow((x - trace(t) * unit_symmetric_tensor<dim>()).norm() <
                  1e-15 * t.norm(),
                ExcInternalError());

    deallog << "x=" << std::endl;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        deallog << i << ' ' << j << ' ' << x[i][j] << std::endl;
  }

  // test 2: check outer product of t with
  // itself
  {
    const SymmetricTensor<4, dim> T = outer_product<dim>(t, t);

    // T*t should yield norm(t)^2*t
    SymmetricTensor<2, dim> x = T * t;
    AssertThrow((x - (t * t) * t).norm() < 1e-15 * t.norm(),
                ExcInternalError());

    deallog << "x=" << std::endl;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        deallog << i << ' ' << j << ' ' << x[i][j] << std::endl;
  }
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
