// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test Utilities::inverse_Hilbert_space_filling_curve in 1D

#include <deal.II/base/utilities.h>

#include "../tests.h"


void
test()
{
  const std::vector<Point<1>> points = {
    Point<1>(1), Point<1>(0), Point<1>(13), Point<1>(-2), Point<1>(22)};

  const int  bit_depth = 5;
  const auto res =
    Utilities::inverse_Hilbert_space_filling_curve(points, bit_depth);

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

  for (const auto ind : index)
    deallog << points[ind][0] << std::endl;
}

int
main()
{
  initlog();

  test();
}
