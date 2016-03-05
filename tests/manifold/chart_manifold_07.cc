//-------------------------------------------------------------------
//    Copyright (C) 2016 by the deal.II authors.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------


// Test direction vector on a cylinder surface

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/manifold.h>


Point<3> periodicity (/*r=*/0,
                            /*phi=*/2*numbers::PI,
                            /*z=*/0);

class MyCylinderManifold : public ChartManifold<2,3,3>
{
public:
  static const int dim = 2;
  static const int spacedim = 3;
  static const int chartdim = 3;

  MyCylinderManifold ()
    :
    ChartManifold<dim,spacedim,spacedim>(periodicity)
  {}


  virtual
  Point<spacedim>
  pull_back(const Point<spacedim> &space_point) const
  {
    const double x = space_point[0];
    const double y = space_point[1];
    const double z = space_point[2];

    const double r   = std::sqrt(x*x + y*y);
    const double phi = std::atan2(y,x);

    return Point<3>(r,
                    phi,
                    z);
  }


  virtual
  Point<spacedim>
  push_forward(const Point<spacedim> &chart_point) const
  {
    const double r   = chart_point[0];
    const double phi = chart_point[1];
    const double z   = chart_point[2];

    return Point<3>(r*std::cos(phi),
                    r*std::sin(phi),
                    z);
  }

  virtual
  DerivativeForm<1,spacedim,spacedim>
  push_forward_gradient(const Point<spacedim> &chart_point) const
  {
    DerivativeForm<1,spacedim,spacedim> g;

    const double r   = chart_point[0];
    const double phi = chart_point[1];
    const double z   = chart_point[2];

    g[0][0] = std::cos(phi);
    g[0][1] = -r*std::sin(phi);
    g[0][2] = 0;

    g[1][0] = std::sin(phi);
    g[1][1] = r*std::cos(phi);
    g[1][2] = 0;

    g[2][0] = 0;
    g[2][1] = 0;
    g[2][2] = 1;

    return g;
  }
};



void test_direction (const Point<3> &x1,
                     const Point<3> &x2)
{
  static MyCylinderManifold manifold;

  // check both the direction x1->x2 and x2->x1
  deallog << '[' << x1 << "] -> [" << x2 << "]: "
          << manifold.get_tangent_vector (x1, x2) << std::endl;
  deallog << '[' << x2 << "] -> [" << x1 << "]: "
          << manifold.get_tangent_vector (x2, x1) << std::endl;
}


void test()
{
  MyCylinderManifold manifold;

  // check two points that are straight up and down from each other
  test_direction (manifold.push_forward (Point<3>(/*r  =*/4,
                                                          /*phi=*/0.1,
                                                          /*z  =*/-1)),
                  manifold.push_forward (Point<3>(/*r  =*/4,
                                                          /*phi=*/0.1,
                                                          /*z  =*/+2)));

  // check two points that are radial
  test_direction (manifold.push_forward (Point<3>(/*r  =*/1,
                                                          /*phi=*/0.1,
                                                          /*z  =*/-1)),
                  manifold.push_forward (Point<3>(/*r  =*/4,
                                                          /*phi=*/0.1,
                                                          /*z  =*/-1)));

  // check two points that are horizontal
  test_direction (manifold.push_forward (Point<3>(/*r  =*/4,
                                                          /*phi=*/0,
                                                          /*z  =*/-1)),
                  manifold.push_forward (Point<3>(/*r  =*/4,
                                                          /*phi=*/numbers::PI/4,
                                                          /*z  =*/-1)));

  // same but rotated
  test_direction (manifold.push_forward (Point<3>(/*r  =*/4,
                                                          /*phi=*/numbers::PI/4,
                                                          /*z  =*/-1)),
                  manifold.push_forward (Point<3>(/*r  =*/4,
                                                          /*phi=*/numbers::PI/2,
                                                          /*z  =*/-1)));


  // check two points that are at the same radius but not horizontal
  test_direction (manifold.push_forward (Point<3>(/*r  =*/4,
                                                          /*phi=*/0,
                                                          /*z  =*/-1)),
                  manifold.push_forward (Point<3>(/*r  =*/4,
                                                          /*phi=*/numbers::PI/4,
                                                          /*z  =*/1)));
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  test();

  return 0;
}

