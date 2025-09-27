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


// compute the deviator tensor as stated in the documentation of outer_product

#include <deal.II/base/symmetric_tensor.h>

#include "../tests.h"


template <int dim>
void
test()
{
  deallog << "dim=" << dim << std::endl;

  const SymmetricTensor<4, dim> T =
    (identity_tensor<dim>() - 1. / dim *
                                outer_product(unit_symmetric_tensor<dim>(),
                                              unit_symmetric_tensor<dim>()));

  AssertThrow((T - deviator_tensor<dim>()).norm() <= 1e-15 * T.norm(),
              ExcInternalError());
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
