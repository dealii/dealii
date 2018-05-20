// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
//---------------------------------------------------------------------

// Create a list of points, and compute the closest points to a given one.

#include "../tests.h"
#include <deal.II/numerics/kdtree.h>

int
main()
{
  initlog();

  KDTree<2> kdtree;

  std::vector<Point<2>> points;

  // Add four points
  points.push_back(Point<2>(0, 0));
  points.push_back(Point<2>(0, 1));
  points.push_back(Point<2>(1, 0));
  points.push_back(Point<2>(1, 1));

  std::vector<Point<2>> test_points;
  test_points.push_back(Point<2>(.5, .5));
  test_points.push_back(Point<2>(2, 0));
  test_points.push_back(Point<2>(2, 2));

  kdtree.set_points(points);

  std::vector<unsigned int> indices;
  std::vector<double>       distances;

  // Get closest points. Do a few rounds
  for(auto& p : test_points)
    {
      for(unsigned int i = 1; i < points.size() + 1; ++i)
        {
          auto res = kdtree.get_closest_points(p, i);

          deallog << std::endl
                  << "The first " << i << " closest points to " << p
                  << " are:" << std::endl;
          for(unsigned int j = 0; j < i; ++j)
            {
              deallog << "points[" << res[j].first
                      << "] = " << points[res[j].first]
                      << ", distance: " << res[j].second << std::endl;
            }
        }
    }
}
