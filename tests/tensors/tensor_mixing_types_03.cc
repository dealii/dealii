// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

// Test for mixed Number type assignments and casts of nested Tensors.

#include <deal.II/base/tensor.h>

#include <complex>

#include "../tests.h"


int
main()
{
  initlog();

  Tensor<1, 1, Tensor<1, 2, float>>  nested_f;
  Tensor<1, 1, Tensor<1, 2, double>> nested_d;

  nested_d = nested_f;
  nested_f = nested_d;

  nested_d = static_cast<Tensor<1, 1, Tensor<1, 2, double>>>(nested_f);
  nested_f = static_cast<Tensor<1, 1, Tensor<1, 2, float>>>(nested_d);

  deallog << "OK!" << std::endl;
}
