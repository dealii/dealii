// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test Tensor<1,dim> * SymmetricTensor<2,dim>

#include <deal.II/base/symmetric_tensor.h>

#include "../tests.h"


template <int dim>
void
test()
{
  deallog << "dim=" << dim << std::endl;

  SymmetricTensor<2, dim> s;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      s[i][j] = (i + 1) * (j + 1);

  Tensor<1, dim> t;
  for (unsigned int i = 0; i < dim; ++i)
    t[i] = (i == 0 ? 1 : 2);

  deallog << s * t << std::endl << t * s << std::endl;

  Assert(s * t == t * s, ExcInternalError());
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test<3>();
  test<5>();

  deallog << "OK" << std::endl;
}
