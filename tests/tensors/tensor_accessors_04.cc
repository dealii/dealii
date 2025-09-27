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

// Test TensorAccessors::contract3

#include <deal.II/base/table_indices.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_accessors.h>

#include "../tests.h"


int
main()
{
  initlog();

  // Contract rank 0, rank 0, rank 0:
  {
    int a = 1;
    int b = 2;
    int c = 4;

    deallog << "left:   " << a << std::endl;
    deallog << "middle: " << b << std::endl;
    deallog << "right:  " << c << std::endl;
    deallog << "Result: " << TensorAccessors::contract3<0, 0, 5, int>(a, b, c)
            << std::endl;
    deallog << std::endl;
  }

  int               c_left[5] = {0, 1, 2, 3, 4};
  Tensor<1, 5, int> left(c_left);

  int                     c_right[5] = {0, 100, 200, 300, 400};
  const Tensor<1, 5, int> right(c_right);


  // Contract rank 1, rank 1, rank 0:
  {
    const int c = -10;

    deallog << "left:   " << left << std::endl;
    deallog << "middle: " << right << std::endl;
    deallog << "right:  " << c << std::endl;
    deallog << "Result: "
            << TensorAccessors::contract3<1, 0, 5, int>(left, right, c)
            << std::endl;
    deallog << std::endl;
  }


  // Contract rank 0, rank 1, rank 1:
  {
    const int c = -5;

    deallog << "left:   " << c << std::endl;
    deallog << "middle: " << left << std::endl;
    deallog << "right:  " << right << std::endl;
    deallog << "Result: "
            << TensorAccessors::contract3<0, 1, 5, int>(c, left, right)
            << std::endl;
    deallog << std::endl;
  }

  Tensor<2, 5, int> middle;
  middle[0] = left;
  middle[1] = 2 * left;
  middle[2] = 4 * left;
  middle[3] = 8 * left;
  middle[4] = 16 * left;

  // Contract rank 1, rank 2, rank 1:
  {
    deallog << "left:   " << left << std::endl;
    deallog << "middle: " << middle << std::endl;
    deallog << "right:  " << right << std::endl;
    deallog << "Result: "
            << TensorAccessors::contract3<1, 1, 5, int>(left, middle, right)
            << std::endl;
    deallog << std::endl;
    // manually verified to be equal to the old implementation
    // deallog << contract3(left, middle, right) << std::endl;
  }

  Tensor<3, 5, int> middle3;
  middle3[0] = middle;
  middle3[1] = 2 * middle;
  middle3[2] = 4 * middle;
  middle3[3] = 8 * middle;
  middle3[4] = 16 * middle;

  // Contract rank 2, rank 3, rank 1:
  {
    deallog << "left:   " << middle << std::endl;
    deallog << "middle: " << middle3 << std::endl;
    deallog << "right:  " << right << std::endl;
    deallog << "Result: "
            << TensorAccessors::contract3<2, 1, 5, int>(middle, middle3, right)
            << std::endl;
    deallog << std::endl;
    //     manually verified to be equal to the old implementation
    //     deallog << contract3(middle, middle3, right) << std::endl;
  }

  // Contract rank 1, rank 3, rank 2:
  {
    deallog << "left:   " << left << std::endl;
    deallog << "middle: " << middle3 << std::endl;
    deallog << "right:  " << middle << std::endl;
    deallog << "Result: "
            << TensorAccessors::contract3<1, 2, 5, int>(left, middle3, middle)
            << std::endl;
    deallog << std::endl;
    //     manually verified to be equal to the old implementation
    //     deallog << contract3(left, middle3, middle) << std::endl;
  }


  Tensor<4, 5, int> middle4;
  middle4[0] = middle3;
  middle4[1] = 2 * middle3;
  middle4[2] = 4 * middle3;
  middle4[3] = 8 * middle3;
  middle4[4] = 16 * middle3;

  Tensor<2, 5, int> left2 = 17 * middle;

  // Contract rank 2, rank 4, rank 2:
  {
    deallog << "left:   " << left2 << std::endl;
    deallog << "middle: " << middle4 << std::endl;
    deallog << "right:  " << middle << std::endl;
    deallog << "Result: "
            << TensorAccessors::contract3<2, 2, 5, long int>(left2,
                                                             middle4,
                                                             middle)
            << std::endl;
    deallog << std::endl;
    //     manually verified to be equal to the old implementation
    //     deallog << contract3(left2, middle4, middle) << std::endl;
  }
}
