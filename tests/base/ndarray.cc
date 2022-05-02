// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


// test that the ndarray type alias does what it promises to do.

#include <deal.II/base/ndarray.h>

#include <type_traits>

#include "../tests.h"

int
main()
{
  initlog();

  dealii::ndarray<double>             r0;
  dealii::ndarray<double, 1>          r1;
  dealii::ndarray<double, 1, 2>       r2;
  dealii::ndarray<double, 1, 2, 3>    r3;
  dealii::ndarray<double, 1, 2, 3, 4> r4;

  static_assert(std::is_same<decltype(r0), double>::value,
                "types must be the same");
  static_assert(std::is_same<decltype(r1), std::array<double, 1>>::value,
                "types must be the same");
  static_assert(
    std::is_same<decltype(r2), std::array<std::array<double, 2>, 1>>::value,
    "types must be the same");
  static_assert(
    std::is_same<decltype(r3),
                 std::array<std::array<std::array<double, 3>, 2>, 1>>::value,
    "types must be the same");
  static_assert(
    std::is_same<decltype(r4),
                 std::array<std::array<std::array<std::array<double, 4>, 3>, 2>,
                            1>>::value,
    "types must be the same");

  return 0;
}
