// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2022 by the deal.II authors
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


// test Utilities::pack_integers() for indices from
// Utilities::inverse_Hilbert_space_filling_curve in 2D on the following points
// with bits=2:

//      3 |    5---6   9---10
//        |    |   |   |   |
//      2 |    4   7---8   11
//        |    |           |
//      1 |    3---2   13--12
//        |        |   |
//      0 |    0---1   14--15
//        |
//         ------------------
//             0   1   2   3

#include <deal.II/base/utilities.h>

#include "../tests.h"


void
test()
{
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
                                        Point<2>(3, 3),
                                        Point<2>(3, 2),
                                        Point<2>(3, 1),
                                        Point<2>(2, 1),
                                        Point<2>(2, 0),
                                        Point<2>(3, 0)};

  const int  bit_depth = 2;
  const auto res =
    Utilities::inverse_Hilbert_space_filling_curve(points, bit_depth);

  for (const auto &p : res)
    deallog << p[0] << ' ' << p[1] << ' '
            << Utilities::pack_integers<2>(p, bit_depth) << std::endl;
}

int
main()
{
  initlog();

  test();
}
