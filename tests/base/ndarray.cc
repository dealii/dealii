// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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

  static_assert(std::is_same_v<decltype(r0), double>, "types must be the same");
  static_assert(std::is_same_v<decltype(r1), std::array<double, 1>>,
                "types must be the same");
  static_assert(
    std::is_same_v<decltype(r2), std::array<std::array<double, 2>, 1>>,
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
