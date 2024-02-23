// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Makes a OpenCASCADE cicrular arc, and project a few points onto it.


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/opencascade/utilities.h>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <GC_MakeCircle.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

#include "../tests.h"

using namespace OpenCASCADE;


int
main()
{
  initlog();

  // A unit circle
  gp_Dir        z_axis(0., 0., 1.);
  gp_Pnt        center(0., 0., 0.);
  gp_Ax2        axis(center, z_axis);
  Standard_Real radius(1.);

  GC_MakeCircle make_circle(axis, radius);
  Handle(Geom_Circle) circle = make_circle.Value();
  TopoDS_Edge edge           = BRepBuilderAPI_MakeEdge(circle);

  // Now get a few points and project
  // them on the circle
  std::vector<Point<3>> points;

  points.push_back(Point<3>(3, 0, 0));
  points.push_back(Point<3>(0, 3, 0));
  // This one is tricky... the point
  // is not on the xy plane. If we
  // put it on the axis (x=0, y=0),
  // there should not exist a
  // projection. Any value for z
  // should give the same result.
  points.push_back(Point<3>(.1, 0, 3));
  points.push_back(Point<3>(.1, 0, 4));


  double       u, v;
  TopoDS_Shape sh;
  for (unsigned int i = 0; i < points.size(); ++i)
    {
      std::tuple<Point<3>, TopoDS_Shape, double, double> ref =
        project_point_and_pull_back(edge, points[i]);

      Point<3> pp = std::get<0>(ref);
      sh          = std::get<1>(ref);
      u           = std::get<2>(ref);
      v           = std::get<3>(ref);

      deallog << "Origin: " << points[i] << ", on unit circle: " << pp
              << ", with local coordinates (u, v): (" << u << ", " << v << ')'
              << std::endl;
    }
  return 0;
}
