// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2009 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test multiplication with a Tensor<1,dim>

#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"


template <int dim>
void
check()
{
  SymmetricTensor<2, dim> S;
  for (unsigned int i = 0; i < S.n_independent_components; ++i)
    S[S.unrolled_to_component_indices(i)] = Testing::rand() % 10;

  Tensor<1, dim> x;
  for (unsigned int i = 0; i < dim; ++i)
    x[i] = Testing::rand() % 10;

  deallog << "S = " << S << std::endl;
  deallog << "x = " << x << std::endl;
  deallog << "S*x = " << S * x << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check<1>();
  check<2>();
  check<3>();

  deallog << "OK" << std::endl;
}
