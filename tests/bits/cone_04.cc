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
// ---------------------------------------------------------------------



// check ConeBoundary<3>::normal_vector()

#include "../tests.h"
#include <deal.II/grid/tria_boundary_lib.h>




void check ()
{
  const Point<3> p1, p2(0,0,2);
  const ConeBoundary<3> boundary (0, 1, p1, p2);

  // test a bunch of points that are randomly chosen. first (manually)
  // project them onto the cone (perpendicular to the z-axis) and then
  // compute the normal vector. we test that they are correct by
  // making sure that they (i) have unit length, (ii) are
  // perpendicular to the tangent vector to the cone that runs from
  // the vertex (at the origin) to the projected point, and (iii) are
  // perpendicular to the tangent vector that is tangent to the circle
  // the projected point sits on
  Point<3> points[] = { Point<3>(1,1,1),
                        Point<3>(1,0,1),
                        Point<3>(0,1,1),

                        Point<3>(1,1,1.5),
                        Point<3>(1,0,1.5),
                        Point<3>(0,1,1.5)
                      };
  for (unsigned int i=0; i<sizeof(points)/sizeof(points[0]); ++i)
    {
      const Point<2> radial_component (points[i][0],
                                       points[i][1]);
      const Point<2> projected_radial_component
        = radial_component/radial_component.norm()
          * (points[i][2]/2);  // radius of cone at given z

      const Point<3> projected_point (projected_radial_component[0],
                                      projected_radial_component[1],
                                      points[i][2]);

      const Tensor<1,3> tangent_1 = projected_point / projected_point.norm();
      const Tensor<1,3> tangent_2 = cross_product_3d (tangent_1,
                                                      Point<3>(0,0,1));

      // get the normal vector and test it
      const Tensor<1,3> normal_vector
        = boundary.normal_vector( /* dummy= */ Triangulation<3>::face_iterator(),
                                               projected_point);

      Assert (std::fabs(normal_vector.norm()-1) < 1e-12,
              ExcInternalError());
      Assert (std::fabs(normal_vector * tangent_1) < 1e-12,
              ExcInternalError());
      Assert (std::fabs(normal_vector * tangent_2) < 1e-12,
              ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


int main ()
{
  initlog();

  check ();
}
