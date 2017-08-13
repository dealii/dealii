//-----------------------------------------------------------
//
//    Copyright (C) 2015 by the deal2lkit authors
//
//    This file is part of the deal2lkit library.
//
//    The deal2lkit library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE at
//    the top level of the deal2lkit distribution.
//
//-----------------------------------------------------------

// Create a list of points, and compute the minimum distance from some other points
// to this set, using a kdtree library

#include "../tests.h"
#include <deal.II/numerics/kdtree_distance.h>


int main ()
{
  initlog();

  KDTreeDistance<2> kdtree;

  std::vector<Point<2> > points;

  // Add four points
  points.push_back(Point<2>(0,0));
  points.push_back(Point<2>(0,1));
  points.push_back(Point<2>(1,0));
  points.push_back(Point<2>(1,1));

  deallog << "Distance from unit square:" << std::endl;

  std::vector<Point<2> > test_points;
  test_points.push_back(Point<2>(.5, .5));
  test_points.push_back(Point<2>(2, 0));
  test_points.push_back(Point<2>(2, 2));

  kdtree.set_points(points);

  for (auto &p : test_points)
    {
      auto res = kdtree.get_closest_points(p,1)[0];
      deallog << "P: " << p << ", distance: " << res.second << ", index: " << res.first << std::endl;
    }

  deallog << "Consistency checking: the following are all the points in the set." << std::endl;
  for (auto &p : points)
    {
      auto res = kdtree.get_closest_points(p,1)[0];
      deallog << "P: " << p << ", distance: " << res.second << ", index: " << res.first << std::endl;
    }
}
