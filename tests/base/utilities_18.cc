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


// test Utilities::inverse_Hilbert_space_filling_curve in (effectively) 1D
// take utilities_16 input and put on a coordinate line in 3D.
// output is made exactly the same as in utilities_16

#include <deal.II/base/utilities.h>

#include "../tests.h"


void
test(unsigned int plane = 1)
{
  Assert(plane < 3, ExcInternalError());
  const std::vector<Point<1>> points = {
    Point<1>(1), Point<1>(0), Point<1>(13), Point<1>(-2), Point<1>(22)};

  std::vector<Point<3>> points_degenerate(points.size());
  for (unsigned int i = 0; i < points.size(); ++i)
    points_degenerate[i][plane] = points[i][0];

  const int  bit_depth = 5;
  const auto res =
    Utilities::inverse_Hilbert_space_filling_curve(points_degenerate,
                                                   bit_depth);

  std::vector<unsigned int> index(points.size());
  for (unsigned int i = 0; i < index.size(); ++i)
    index[i] = i;

  std::sort(index.begin(),
            index.end(),
            [&](const unsigned int &a, const unsigned int &b) {
              AssertIndexRange(a, res.size());
              AssertIndexRange(b, res.size());
              return std::lexicographical_compare(res[a].begin(),
                                                  res[a].end(),
                                                  res[b].begin(),
                                                  res[b].end());
            });

  deallog << plane << ':' << std::endl;
  for (const auto ind : index)
    {
      AssertThrow(res[ind][0] == 0, ExcInternalError());
      AssertThrow(res[ind][1] == 0, ExcInternalError());
      deallog << points[ind][0] << std::endl;
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
