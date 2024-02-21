// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// check that operations on tensors of vectorized arrays are properly
// supported

#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

#include "../tests.h"

int
main()
{
  initlog();

  Tensor<1, 2, VectorizedArray<double>> t;
  deallog << "Tensor norm: " << t.norm()[0] << std::endl;
  Point<3, VectorizedArray<float>> p;
  deallog << "Point norm: " << p.norm()[0] << std::endl;
  SymmetricTensor<2, 3, VectorizedArray<double>> st2;
  SymmetricTensor<4, 2, VectorizedArray<double>> st4;
  deallog << "Symmetric tensor norm: " << st2.norm()[0] << ' ' << st4.norm()[0]
          << std::endl;
}
