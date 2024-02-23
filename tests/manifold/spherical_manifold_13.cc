// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test that the algorithm used in SphericalManifold<3>::get_new_point is stable
// with respect to distortion of the points and weights.

#include <deal.II/base/point.h>

#include <deal.II/grid/manifold_lib.h>

#include <random>

#include "../tests.h"

int
main()
{
  initlog();

  const int                       n_iterations = 1000;
  std::mt19937                    generator;
  std::uniform_int_distribution<> distribution(-1., 1.);

  SphericalManifold<3> manifold;

  for (int j = 2; j < 20; ++j)
    {
      const double eps = std::pow(10., -j);

      std::vector<Point<3>> points_ref{{-.1, -.1, .5},
                                       {.1, -.1, .5},
                                       {-.1, .1, .5},
                                       {.1, .1, .5}};
      std::vector<double>   weights_ref{.25, .25, .25, .25};

      const Point<3> new_point_ref =
        manifold.get_new_point(make_array_view(points_ref),

                               make_array_view(weights_ref));

      double max_difference = 0.;

      for (unsigned int i = 0; i < n_iterations; ++i)
        {
          std::vector<Point<3>> points  = {{-.1, -.1, .5},
                                           {.1, -.1, .5},
                                           {-.1, .1, .5},
                                           {.1, .1, .5}};
          std::vector<double>   weights = {.25, .25, .25, .25};

          for (Point<3> &point : points)
            {
              point[0] += distribution(generator) * eps;
              point[1] += distribution(generator) * eps;
              point[2] += distribution(generator) * eps;
            }

          for (double &weight : weights)
            weight += distribution(generator) * eps;

          const Point<3> new_point =
            manifold.get_new_point(make_array_view(points),
                                   make_array_view(weights));
          const Tensor<1, 3> difference      = new_point - new_point_ref;
          const double       difference_norm = difference.norm();
          max_difference = std::max(difference_norm, max_difference);
        }

      if (max_difference / eps > 5)
        deallog << "Distortion: " << eps
                << " max difference: " << max_difference << std::endl;
      else
        deallog << "Distortion " << eps << " OK" << std::endl;
    }
  return 0;
}
