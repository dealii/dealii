// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test Utilities::inverse_Hilbert_space_filling_curve in 2D
// on the following nine points with bits=4:

//        |
//     15 |    @---@   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |    |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
//        |    @   @---@   @   @   @---@   @   @   @---@   @   @   @---@   @
//        |    |           |   |           |   |           |   |           |
//        |    @---@   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |        |   |           |   |           |   |           |   |
//        |    @---@   @---@---@---@   @---@   @---@   @---@---@---@   @---@
//        |    |                           |   |                           |
//        |    @   @---@---@   @---@---@   @   @   @---@---@   @---@---@   @
//        |    |   |       |   |       |   |   |   |       |   |       |   |
// Axes[1]|    @---@   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |            |           |                   |           |
//        |    @---@   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |    |   |       |   |       |   |   |   |       |   |       |   |
//        |    @   @---@---@   @---@---@   @---@   @---@---@   @---@---@   @
//        |    |                                                           |
//        |    @---@   @---@---@   @---@---@   @---@---@   @---@---@   @---@
//        |        |   |       |   |       |   |       |   |       |   |
//        |    @---@   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |    |           |           |           |           |           |
//        |    @   @---@   @   @---@   @---@   @---@   @---@   @   @---@   @
//        |    |   |   |   |   |   |       |   |       |   |   |   |   |   |
//        |    @---@   @---@   @   @---@---@   @---@---@   @   @---@   @---@
//        |                    |                           |
//      3 |    5---6   9---@   @   @---@---@   @---@---@   @   @---@   @---@
//        |    |   |   |   |   |   |       |   |       |   |   |   |   |   |
//      2 |    4   7---8   @   @---@   @---@   @---@   @---@   @   @---@   @
//        |    |           |           |           |           |           |
//      1 |    3---2   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |        |   |       |   |       |   |       |   |       |   |
//      0 |    0---1   @---@---@   @---@---@   @---@---@   @---@---@   @--255
//        |
//         -------------------------------------------------------------------
//             0   1   2   3          ---> Axes[0]                         15


#include <deal.II/base/utilities.h>

#include "../tests.h"


void
test1()
{
  deallog << "test1:" << std::endl;
  const std::vector<Point<2>> points = {Point<2>(0, 0),
                                        Point<2>(1, 0),
                                        Point<2>(1, 1),
                                        Point<2>(0, 1),
                                        Point<2>(0, 2),
                                        Point<2>(0, 3),
                                        Point<2>(1, 3),
                                        Point<2>(1, 2),
                                        Point<2>(2, 2),
                                        Point<2>(2, 3),
                                        Point<2>(15, 15),
                                        Point<2>(15, 0)};

  const auto res = Utilities::inverse_Hilbert_space_filling_curve(points, 4);

  for (const auto &p : res)
    deallog << p[0] << ' ' << p[1] << std::endl;
}

// https://github.com/aditi137/Hilbert/blob/master/Hilbert/test.py
// test int version in 3D
void
test2()
{
  deallog << "test2:" << std::endl;
  std::vector<std::array<std::uint64_t, 3>> points = {{{5, 10, 20}}};

  const auto res = Utilities::inverse_Hilbert_space_filling_curve<3>(points, 5);

  for (const auto &p : res)
    deallog << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
}


// finally, a quick test of fractional coordinates -> int conversion inside
// inverse_Hilbert_space_filling_curve():
template <typename Integer, typename Number>
void
test3()
{
  const int bits_int = std::numeric_limits<Integer>::digits;
  const int bits_num = std::numeric_limits<Number>::digits;

  const int min_bits = std::min(bits_int, bits_num);

  const Integer max =
    (min_bits == bits_int ? std::numeric_limits<Integer>::max() :
                            ((Integer)1 << min_bits) - 1);

  // make sure conversion is within bounds:
  {
    const Number  v = 0.2;
    const Integer i = (Integer)(v * (Number)max);
    Assert(i > 0 && i < max, ExcInternalError());
  }

  {
    const Number  v = 0.7;
    const Integer i = (Integer)(v * (Number)max);
    Assert(i > 0 && i < max, ExcInternalError());
  }

  {
    const Number  v = 0.;
    const Integer i = (Integer)(v * (Number)max);
    Assert(i == 0, ExcInternalError());
  }

  {
    const Number  v = 1.;
    const Integer i = (Integer)(v * (Number)max);
    Assert(i == max, ExcInternalError());
  }
}

int
main()
{
  initlog();

  test1();
  test2();

  test3<std::uint64_t, long double>();
}
