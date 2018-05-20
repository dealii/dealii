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

// Create a list of points, and the ones within a specified radius of the target points

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
  test_points.push_back(Point<2>(2, 0));
  test_points.push_back(Point<2>(2, 2));

  std::vector<double> radii;
  radii.push_back(.8);
  radii.push_back(1.1);
  radii.push_back(1.5);
  radii.push_back(2);

  kdtree.set_points(points);

  std::vector<std::pair<unsigned int, double>> matches;

  // Get points within ball
  for(auto& p : test_points)
    {
      for(auto& r : radii)
        {
          auto matches = kdtree.get_points_within_ball(p, r, true);

          deallog << std::endl
                  << "At distance less than " << r << " from " << p
                  << " we have:" << std::endl;
          for(auto& m : matches)
            {
              deallog << "points[" << m.first << "] = " << points[m.first]
                      << ", distance: " << m.second << std::endl;
            }
        }
    }
}
