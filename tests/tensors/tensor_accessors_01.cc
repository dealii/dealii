// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_accessors.h>

#include "../tests.h"

#define PRINTME(bar)                       \
  for (unsigned int i = 0; i < 2; ++i)     \
    for (unsigned int j = 0; j < 2; ++j)   \
      for (unsigned int k = 0; k < 2; ++k) \
        deallog << bar[i][j][k] << ' ';    \
  deallog << std::endl;


int
main()
{
  initlog();

  Tensor<9, 3, int> t;
  t[0][1][2][0][1][2][0][1][2] = 42;

  // Reorder index 4 (count begins at 0) to last place:
  TensorAccessors::internal::
    ReorderedIndexView<4, 9, Tensor<9, 3, int>> // auto ...
      foo = TensorAccessors::reordered_index_view<4, 9>(t);

  // test access and assignment:
  {
    //               0  1  2  3  5  6  7  8  4
    deallog << foo[0][1][2][0][2][0][1][2][1] << std::endl;

    int  temp = foo[0][1][2][0][2][0][1][2][1];
    int &ref  = foo[0][1][2][0][2][0][1][2][1];
    ref       = temp;

    foo[0][1][2][0][2][0][1][2][1] = temp + ref;

    deallog << t[0][1][2][0][1][2][0][1][2] << std::endl;
  }

  // test read-only access:
  {
    const Tensor<9, 3, int> &t_ref = t;

    TensorAccessors::internal::
      ReorderedIndexView<4, 9, const Tensor<9, 3, int>> // auto ...
        const_foo = TensorAccessors::reordered_index_view<4, 9>(t_ref);

    //                     0  1  2  3  5  6  7  8  4
    deallog << const_foo[0][1][2][0][2][0][1][2][1] << std::endl;

    int &tmp = foo[0][1][2][0][2][0][1][2][1];
    tmp      = const_foo[0][1][2][0][2][0][1][2][1] / 2;
    deallog << const_foo[0][1][2][0][2][0][1][2][1] << std::endl;
  }

  // test nested reordering:
  {
    TensorAccessors::internal::ReorderedIndexView<
      0,
      9,
      TensorAccessors::internal::
        ReorderedIndexView<4, 9, Tensor<9, 3, int>>> // auto ...
      foo2 = TensorAccessors::reordered_index_view<0, 9>(foo);

    //              t 0  1  2  3  4  5  6  7  8
    //         access 0  1  2  0  1  2  0  1  2

    //           foo  0  1  2  3  5  6  7  8  4
    //           foo2 1  2  3  5  6  7  8  4  0
    deallog << foo2[1][2][0][2][0][1][2][1][0] << std::endl;
  }

  {
    // check whether all special cases of reordering work as expected:

    int initializer[2][2][2] = {{{0, 1}, {2, 3}}, {{4, 5}, {6, 7}}};
    Tensor<3, 2, int> t(initializer);

    deallog << "Order of indices 0 1 2  -->  ";
    TensorAccessors::internal::
      ReorderedIndexView<2, 3, Tensor<3, 2, int>> // auto ...
        foo012 = TensorAccessors::reordered_index_view<2, 3>(t);
    PRINTME(foo012);

    deallog << "Order of indices 0 2 1  -->  ";
    TensorAccessors::internal::
      ReorderedIndexView<1, 3, Tensor<3, 2, int>> // auto ...
        foo021 = TensorAccessors::reordered_index_view<1, 3>(t);
    PRINTME(foo021);

    deallog << "Order of indices 1 2 0  -->  ";
    TensorAccessors::internal::
      ReorderedIndexView<0, 3, Tensor<3, 2, int>> // auto ...
        foo120 = TensorAccessors::reordered_index_view<0, 3>(t);
    PRINTME(foo120);

    deallog << "Order of indices 1 0 2  -->  ";
    TensorAccessors::internal::ReorderedIndexView<
      1,
      3,
      TensorAccessors::internal::
        ReorderedIndexView<0, 3, Tensor<3, 2, int>>> // auto ...
      foo102 = TensorAccessors::reordered_index_view<1, 3>(foo120);
    PRINTME(foo102);

    deallog << "Order of indices 2 0 1  -->  ";
    TensorAccessors::internal::ReorderedIndexView<
      0,
      3,
      TensorAccessors::internal::
        ReorderedIndexView<0, 3, Tensor<3, 2, int>>> // auto ...
      foo201 = TensorAccessors::reordered_index_view<0, 3>(foo120);
    PRINTME(foo201);

    deallog << "Order of indices 2 1 0  -->  ";
    TensorAccessors::internal::ReorderedIndexView<
      0,
      3,
      TensorAccessors::internal::
        ReorderedIndexView<1, 3, Tensor<3, 2, int>>> // auto ...
      foo210 = TensorAccessors::reordered_index_view<0, 3>(foo021);
    PRINTME(foo210);
  }

  {
    // check with c-style arrays:
    double t[3][3][3][3][3];
    t[0][1][2][0][1] = 42.;

    dealii::TensorAccessors::internal::
      ReorderedIndexView<2, 5, double[3][3][3][3][3]> // auto ...
        foo = TensorAccessors::reordered_index_view<2, 5>(t);
    deallog << foo[0][1][0][1][2] << std::endl;

    const double(&t_ref)[3][3][3][3][3] = t;
    dealii::TensorAccessors::internal::
      ReorderedIndexView<2, 5, const double[3][3][3][3][3]> // auto ...
        foo2 = TensorAccessors::reordered_index_view<2, 5>(t_ref);
    deallog << foo2[0][1][0][1][2] << std::endl;
  }

  {
    // Is it possible to call the simplest case (where we have to do
    // absolutely nothing?

    Tensor<1, 3, int> t;
    TensorAccessors::internal::
      ReorderedIndexView<0, 1, Tensor<1, 3, int>> // auto ...
        foo = TensorAccessors::reordered_index_view<0, 1>(t);
  }
}
