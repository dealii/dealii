// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// check Tensor::operator= (double)

#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  double a[3][3] = {{1, 2, 3}, {3, 4, 5}, {6, 7, 8}};

  const unsigned int dim = 3;
  Tensor<2, dim>     t(a);

  deallog << t.norm() << std::endl;
  t = 0;
  deallog << t.norm() << std::endl;

  Assert(t.norm() == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}
