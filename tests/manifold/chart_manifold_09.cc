// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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
// ---------------------------------------------------------------------

// Check ChartManifold for periodicity issues. Check that the it finds
// the right intermediate points independent on the number of
// surrounding points

#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/grid/manifold_lib.h>


int
main()
{
  initlog();

  const FunctionManifold<2,2,2> manifold("x;y", "x;y", Point<2>(1.0, 0.0));

  // Some points that would cross the periodicity boundary
  std::vector<Point<2> > points;

  points.push_back(Point<2>(0.1, 0.0));
  points.push_back(Point<2>(0.9, 0.0));

  // And the weights
  std::vector<double> weights(2);

  unsigned int n_intermediates=10;

  // First, use only two surrounding points. The second time around,
  // instead, use four, with additional zero weights
  for (unsigned int test_no=0; test_no<2; ++test_no)
    {
      if (test_no == 1)
        {
          // Test that we get the same result with four surrounding points
          points.push_back(Point<2>(0.3, 1.0));
          points.push_back(Point<2>(0.7, 1.0));
          weights.push_back(0.0);
          weights.push_back(0.0);
        }

      deallog << "Test " << test_no << std::endl;
      for (unsigned int i=0; i<n_intermediates; ++i)
        {
          weights[0] = (double)i/((double)n_intermediates-1);
          weights[1] = 1.0-weights[0];

          deallog << manifold.get_new_point(Quadrature<2>(points, weights)) << std::endl;
        }
    }

  return 0;
}



