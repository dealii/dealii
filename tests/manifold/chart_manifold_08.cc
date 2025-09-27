// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test direction vector on a cylinder surface. this is like the _07
// test, but choose points in such a way that we wrap around the
// periodicity in phi

#include <deal.II/grid/manifold.h>

#include "../tests.h"


Tensor<1, 3>
  periodicity(static_cast<Tensor<1, 3>>(Point<3>(0, 2 * numbers::PI, 0)));

class MyCylinderManifold : public ChartManifold<2, 3, 3>
{
public:
  static const int dim      = 2;
  static const int spacedim = 3;
  static const int chartdim = 3;

  MyCylinderManifold()
    : ChartManifold<dim, spacedim, spacedim>(periodicity)
  {}

  virtual std::unique_ptr<Manifold<2, 3>>
  clone() const override
  {
    return std::unique_ptr<Manifold<2, 3>>(new MyCylinderManifold());
  }

  virtual Point<spacedim>
  pull_back(const Point<spacedim> &space_point) const override
  {
    const double x = space_point[0];
    const double y = space_point[1];
    const double z = space_point[2];

    const double r   = std::sqrt(x * x + y * y);
    const double phi = std::atan2(y, x);

    return Point<3>(r, phi, z);
  }


  virtual Point<spacedim>
  push_forward(const Point<spacedim> &chart_point) const override
  {
    const double r   = chart_point[0];
    const double phi = chart_point[1];
    const double z   = chart_point[2];

    return Point<3>(r * std::cos(phi), r * std::sin(phi), z);
  }

  virtual DerivativeForm<1, spacedim, spacedim>
  push_forward_gradient(const Point<spacedim> &chart_point) const override
  {
    DerivativeForm<1, spacedim, spacedim> g;

    const double r   = chart_point[0];
    const double phi = chart_point[1];
    const double z   = chart_point[2];

    g[0][0] = std::cos(phi);
    g[0][1] = -r * std::sin(phi);
    g[0][2] = 0;

    g[1][0] = std::sin(phi);
    g[1][1] = r * std::cos(phi);
    g[1][2] = 0;

    g[2][0] = 0;
    g[2][1] = 0;
    g[2][2] = 1;

    return g;
  }
};



void
test_direction(const Point<3> &x1, const Point<3> &x2)
{
  static MyCylinderManifold manifold;

  // check both the direction x1->x2 and x2->x1
  deallog << '[' << x1 << "] -> [" << x2
          << "]: " << manifold.get_tangent_vector(x1, x2) << std::endl;
  deallog << '[' << x2 << "] -> [" << x1
          << "]: " << manifold.get_tangent_vector(x2, x1) << std::endl;
}


void
test()
{
  MyCylinderManifold manifold;

  // check two points that are horizontal
  test_direction(manifold.push_forward(Point<3>(/*r  =*/2,
                                                /*phi=*/3 * numbers::PI / 4,
                                                /*z  =*/-1)),
                 manifold.push_forward(Point<3>(/*r  =*/2,
                                                /*phi=*/-3 * numbers::PI / 4,
                                                /*z  =*/-1)));

  // same but rotated
  test_direction(manifold.push_forward(Point<3>(/*r  =*/2,
                                                /*phi=*/-numbers::PI / 4,
                                                /*z  =*/-1)),
                 manifold.push_forward(Point<3>(/*r  =*/2,
                                                /*phi=*/numbers::PI / 4,
                                                /*z  =*/-1)));
}

int
main()
{
  initlog();

  test();

  return 0;
}
