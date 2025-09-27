// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check SphericalManifold for get_intermediate_point and get_tangent_vector
// issues.

#include <deal.II/base/quadrature.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"

int
main()
{
  initlog();

  // // Center and radius of the Ball
  // double radius = center.norm();

  {
    Point<2>                      center(0.0, 0.0);
    const SphericalManifold<2, 2> manifold(center);

    Point<2> P1(1.0, 0.0);
    Point<2> P2(0.0, 1.0);

    Point<2> Q = manifold.get_intermediate_point(P1, P2, .5);

    deallog << "=================================" << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .125) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .25) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .375) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .5) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .625) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .75) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .875) << std::endl;
    deallog << "=================================" << std::endl;
    deallog << manifold.get_intermediate_point(P1, Q, .25) << std::endl;
    deallog << manifold.get_intermediate_point(P1, Q, .5) << std::endl;
    deallog << manifold.get_intermediate_point(P1, Q, .75) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .5) << std::endl;
    deallog << manifold.get_intermediate_point(Q, P2, .25) << std::endl;
    deallog << manifold.get_intermediate_point(Q, P2, .5) << std::endl;
    deallog << manifold.get_intermediate_point(Q, P2, .75) << std::endl;
    deallog << "=================================" << std::endl;
  }

  {
    Point<2>                      center(0.0, 0.0);
    const SphericalManifold<1, 2> manifold(center);

    Point<2> P1(1.0, 0.0);
    Point<2> P2(0.0, 1.0);

    Point<2> Q = manifold.get_intermediate_point(P1, P2, .5);

    deallog << "=================================" << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .125) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .25) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .375) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .5) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .625) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .75) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .875) << std::endl;
    deallog << "=================================" << std::endl;
    deallog << manifold.get_intermediate_point(P1, Q, .25) << std::endl;
    deallog << manifold.get_intermediate_point(P1, Q, .5) << std::endl;
    deallog << manifold.get_intermediate_point(P1, Q, .75) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .5) << std::endl;
    deallog << manifold.get_intermediate_point(Q, P2, .25) << std::endl;
    deallog << manifold.get_intermediate_point(Q, P2, .5) << std::endl;
    deallog << manifold.get_intermediate_point(Q, P2, .75) << std::endl;
    deallog << "=================================" << std::endl;
  }

  {
    Point<3>                      center(0.0, 0.0, 0.0);
    const SphericalManifold<2, 3> manifold(center);

    Point<3> P1(1.0, 0.0, 0.0);
    Point<3> P2(0.0, 0.0, 1.0);

    Point<3> Q = manifold.get_intermediate_point(P1, P2, .5);

    deallog << "=================================" << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .125) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .25) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .375) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .5) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .625) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .75) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .875) << std::endl;
    deallog << "=================================" << std::endl;
    deallog << manifold.get_intermediate_point(P1, Q, .25) << std::endl;
    deallog << manifold.get_intermediate_point(P1, Q, .5) << std::endl;
    deallog << manifold.get_intermediate_point(P1, Q, .75) << std::endl;
    deallog << manifold.get_intermediate_point(P1, P2, .5) << std::endl;
    deallog << manifold.get_intermediate_point(Q, P2, .25) << std::endl;
    deallog << manifold.get_intermediate_point(Q, P2, .5) << std::endl;
    deallog << manifold.get_intermediate_point(Q, P2, .75) << std::endl;
    deallog << "=================================" << std::endl;
  }

  {
    Point<3>                      center(0.0, 0.0, 0.0);
    const SphericalManifold<3, 3> manifold(center);

    Point<3> P1(2.0, 0.0, 0.0);
    Point<3> P2(0.0, std::sqrt(2), std::sqrt(2));

    Point<3> Q = manifold.get_intermediate_point(P1, P2, .5);

    const unsigned int num_points = 20;
    deallog << "=================================" << std::endl;
    for (unsigned int i = 0; i < num_points; ++i)
      deallog << manifold.get_intermediate_point(P1,
                                                 P2,
                                                 (1.0 * i) / (num_points - 1))
              << std::endl;
    deallog << "=================================" << std::endl;
  }

  {
    Point<3>                      center(0.0, 0.0, 0.0);
    const SphericalManifold<3, 3> manifold(center);

    Point<3> P1(1.0, 0.0, 0.0);
    Point<3> P2(0.0, 1.0, 0.0);
    Point<3> P3(0.0, 0.0, 1.0);

    std::vector<Point<3>> points1(3);
    std::vector<Point<3>> points2(3);
    std::vector<Point<3>> points3(3);
    std::vector<double>   weights(3);

    points1[0] = P1;
    points1[1] = P2;
    points1[2] = P3;

    points2[0] = P2;
    points2[1] = P1;
    points2[2] = P3;

    points3[0] = P2;
    points3[1] = P3;
    points3[2] = P1;

    weights[0] = 1.0 / 3.0;
    weights[1] = 1.0 / 3.0;
    weights[2] = 1.0 / 3.0;

    Point<3> Q = manifold.get_new_point(make_array_view(points1),
                                        make_array_view(weights));
    Point<3> S = manifold.get_new_point(make_array_view(points2),
                                        make_array_view(weights));
    Point<3> T = manifold.get_new_point(make_array_view(points3),
                                        make_array_view(weights));

    Point<3> P5(0.707107, 0.707107, 0.0);
    Point<3> P4(0.0, 0.0, 1.0);
    Point<3> R = manifold.get_intermediate_point(P5, P4, 2.0 / 3.0);

    deallog << "=================================" << std::endl;
    deallog << Q << std::endl;
    deallog << S << std::endl;
    deallog << T << std::endl;
    deallog << R << std::endl;
    deallog << "=================================" << std::endl;
  }

  // Quadrature (const std::vector< Point< dim > > &points, const std::vector<
  // double > &weights);
  return 0;
}
