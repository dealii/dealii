// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
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

#include "../tests.h"
#include <deal.II/base/table_indices.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_accessors.h>

#include <array>

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

  const Tensor<9, 3, int>& t_ref = t;

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
