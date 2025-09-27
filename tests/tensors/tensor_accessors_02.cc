// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/table_indices.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_accessors.h>

#include <array>

#include "../tests.h"

int
main()
{
  initlog();

  Tensor<9, 3, int> t;
  t[0][1][2][0][1][2][0][1][2] = 42;

  TableIndices<9> indices(0, 1, 2, 0, 1, 2, 0, 1, 2);

  deallog << TensorAccessors::extract<9>(t, indices) << std::endl;

  TensorAccessors::extract<9>(t, indices) = 84;

  deallog << t[0][1][2][0][1][2][0][1][2] << std::endl;

  const Tensor<9, 3, int> &t_ref = t;

  deallog << TensorAccessors::extract<9>(t_ref, indices) << std::endl;

  // Reduced access:
  {
    deallog << TensorAccessors::extract<8>(t_ref, indices)[2] << std::endl;
    deallog << TensorAccessors::extract<7>(t_ref, indices)[1][2] << std::endl;
    deallog << TensorAccessors::extract<6>(t_ref, indices)[0][1][2]
            << std::endl;
    deallog << TensorAccessors::extract<5>(t_ref, indices)[2][0][1][2]
            << std::endl;
    deallog << TensorAccessors::extract<4>(t_ref, indices)[1][2][0][1][2]
            << std::endl;
    deallog << TensorAccessors::extract<3>(t_ref, indices)[0][1][2][0][1][2]
            << std::endl;
    deallog << TensorAccessors::extract<2>(t_ref, indices)[2][0][1][2][0][1][2]
            << std::endl;
    deallog << TensorAccessors::extract<1>(t_ref,
                                           indices)[1][2][0][1][2][0][1][2]
            << std::endl;
    deallog << TensorAccessors::extract<0>(t_ref,
                                           indices)[0][1][2][0][1][2][0][1][2]
            << std::endl;
  }

  // Apply on c-style arrays:
  {
    double foo[3][3][3][3][3];

    foo[2][1][0][2][1] = 42.;

    unsigned int indices[5];
    indices[0] = 2;
    indices[1] = 1;
    indices[2] = 0;
    indices[3] = 2;
    indices[4] = 1;

    deallog << TensorAccessors::extract<5>(foo, indices) << std::endl;

    // read-only:

    const double(&foo2)[3][3][3][3][3] = foo;
    deallog << TensorAccessors::extract<5>(foo2, indices) << std::endl;

    // via std::array

    // via std::array:
    std::array<unsigned int, 5> temp{{2, 1, 0, 2, 1}};
    deallog << TensorAccessors::extract<5>(foo, temp) << std::endl;
  }
}
