// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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


// Test that a const SymmetricTensor is not assignable, i.e., the accessors
// return a const reference

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

#include "../tests.h"

int
main()
{
  initlog();

  {
    const SymmetricTensor<2, 3, double> st{};
    static_assert(std::is_same_v<decltype(st[0][0]), const double &>);
  }

  {
    const SymmetricTensor<2, 3, VectorizedArray<double>> stv{};
    static_assert(std::is_same<decltype(stv[0][0]),
                               const VectorizedArray<double> &>::value);
  }

  deallog << "OK" << std::endl;
}
