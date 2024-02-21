// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
