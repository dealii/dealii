// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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


// similar to utilities_14, test
// Utilities::inverse_Hilbert_space_filling_curve in 3D for point on the planes
// with bits=2.
// the test output is made to be exactly the same as utilities_14

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
test(const unsigned int plane = 1)
{
  Assert(plane < 3, ExcInternalError());
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

  std::vector<Point<3>> points_degenerate(points.size());
  for (unsigned int i = 0; i < points.size(); ++i)
    {
      unsigned int eff_d = 0;
      for (unsigned int d = 0; d < 3; ++d)
        if (d != plane)
          {
            points_degenerate[i][d] = points[i][eff_d];
            ++eff_d;
          }
    }


  const int  bit_depth = 2;
  const auto res =
    Utilities::inverse_Hilbert_space_filling_curve(points_degenerate,
                                                   bit_depth);

  deallog << plane << ':' << std::endl;
  for (const auto &p : res)
    {
      AssertThrow(p[0] == 0, ExcInternalError());
      deallog << p[1] << ' ' << p[2] << ' '
              << Utilities::pack_integers<3>(p, bit_depth) << std::endl;
    }
  deallog << std::endl;
}

int
main()
{
  initlog();

  test(0);
  test(1);
  test(2);
}
