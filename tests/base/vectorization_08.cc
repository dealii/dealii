// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// check that operations on tensors of vectorized arrays are properly
// supported

#include "../tests.h"

#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

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
  deallog << "Symmetric tensor norm: " << st2.norm()[0] << " " << st4.norm()[0]
          << std::endl;
}
