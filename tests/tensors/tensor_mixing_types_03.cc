// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
