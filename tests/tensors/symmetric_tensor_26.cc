// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


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
