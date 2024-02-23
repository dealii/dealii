// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test TensorAccessors::contract

#include <deal.II/base/table_indices.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_accessors.h>

#include "../tests.h"


int
main()
{
  initlog();

  int               c_left1[5] = {0, 1, 2, 3, 4};
  Tensor<1, 5, int> left1(c_left1);
  int               c_right1[5] = {0, 10, 20, 30, 40};
  Tensor<1, 5, int> right1(c_right1);

  deallog << "left1: " << left1 << std::endl;
  deallog << "right1: " << right1 << std::endl;

  // Apply contract with no_contr == 0, this is a plain outer product:
  {
    Tensor<2, 5, int> result;
    TensorAccessors::contract<0, 1, 1, 5>(result, left1, right1);

    deallog << "Outer product of left1 and right1:" << std::endl;
    deallog << result << std::endl;
  }

  // Apply contract with no_contr == 1, this is a scalar product:
  {
    int result = 0;
    TensorAccessors::contract<1, 1, 1, 5>(result, left1, right1);

    deallog << "Scalar product of left1 and right1:" << std::endl;
    deallog << result << std::endl;
  }

  Tensor<2, 5, int> left2;
  left2[0] = left1;
  left2[1] = 2 * left1;
  left2[2] = 4 * left1;
  left2[3] = 8 * left1;
  left2[4] = 16 * left1;

  Tensor<2, 5, int> right2;
  right2[0] = right1;
  right2[1] = 2 * right1;
  right2[2] = 4 * right1;
  right2[3] = 8 * right1;
  right2[4] = 16 * right1;

  deallog << "left2: " << left2 << std::endl;
  deallog << "right2: " << right2 << std::endl;

  // Apply contract with no_contr == 0, this is a plain outer product:
  {
    Tensor<4, 5, int> result;
    TensorAccessors::contract<0, 2, 2, 5>(result, left2, right2);

    deallog << "Outer product of left2 and right2:" << std::endl;
    deallog << result << std::endl;
    // The result is verified by hand:
    // for (unsigned int i = 0; i < 4; ++i)
    //   for (unsigned int j = 0; j < 4; ++j)
    //     for (unsigned int k = 0; k < 4; ++k)
    //       for (unsigned int l = 0; l < 4; ++l)
    //         {
    //           deallog << '(' << i << ',' << j << ',' << k << ',' << l << "):
    //           "; deallog << result[i][j][k][l] << " = " << left2[i][j] << " *
    //           " << right2[k][l] << std::endl;
    //         }
  }

  // Apply contract with no_contr == 1, and switch indices of right2. This
  // corresponds to a contraction of the last index of left2 with the first
  // index of right2:
  {
    dealii::TensorAccessors::internal::
      ReorderedIndexView<0, 2, dealii::Tensor<2, 5, int>>
        reordered = TensorAccessors::reordered_index_view<0, 2>(right2);

    Tensor<2, 5, int> result;
    TensorAccessors::contract<1, 2, 2, 5>(result, left2, reordered);

    deallog
      << "Contract the last index of left2 with the first index of right2:"
      << std::endl;
    deallog << result << std::endl;
    // Verified to be the same as the old implementation
    // deallog << contract(left2, right2) << std::endl;
  }

  // Apply contract with no_contr == 2, this is a double contraction.
  {
    int result;
    TensorAccessors::contract<2, 2, 2, 5>(result, left2, right2);

    deallog << "Double contraction:" << std::endl;
    deallog << result << std::endl;
    // Verified to be the same as the old implementation
    // deallog << double_contract(left2, right2) << std::endl;
  }
}
