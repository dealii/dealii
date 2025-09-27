// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Create a BSpline surface, and test NormalToMeshProjectionManifold  projection
// using 5 surrounding points.

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>

#include <BRepFill.hxx>
#include <Standard_Stream.hxx>
#include <TopTools.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

#include "../tests.h"

using namespace OpenCASCADE;

int
main()
{
  initlog();

  // Create a bspline passing through the points
  std::vector<Point<3>> pts;
  pts.push_back(Point<3>(0, 0, 0));
  pts.push_back(Point<3>(1, 0, 0));
  // pts.push_back(Point<3>(0.5, 2, -0.5));
  pts.push_back(Point<3>(1, 0, -1));
  pts.push_back(Point<3>(0, 0, -1));

  TopoDS_Edge edge1 = interpolation_curve(pts);
  for (unsigned int i = 0; i < pts.size(); ++i)
    pts[i] += Point<3>(0, 1, 0);
  TopoDS_Edge edge2 = interpolation_curve(pts);

  TopoDS_Face face = BRepFill::Face(edge1, edge2);

  NormalToMeshProjectionManifold<2, 3> manifold(face);
  FlatManifold<2, 3>                   flat_manifold;
  std::vector<Point<3>>                surrounding_points;
  surrounding_points.push_back(Point<3>(0., 0., 0.));
  surrounding_points.push_back(Point<3>(1., 0., 0.));
  Vector<double> weights(surrounding_points.size());
  weights = 1. / surrounding_points.size();
  auto intermediate_point =
    manifold.get_new_point(ArrayView<Point<3>>(&surrounding_points[0],
                                               surrounding_points.size()),
                           ArrayView<double>(&weights[0], weights.size()));
  surrounding_points.push_back(intermediate_point);
  surrounding_points.push_back(Point<3>(0., 1., 0.));
  surrounding_points.push_back(Point<3>(1., 1., 0.));
  weights.reinit(surrounding_points.size());
  weights = 1. / surrounding_points.size();


  auto new_point =
    manifold.get_new_point(ArrayView<Point<3>>(&surrounding_points[0],
                                               surrounding_points.size()),
                           ArrayView<double>(&weights[0], weights.size()));
  auto new_point_flat =
    flat_manifold.get_new_point(ArrayView<Point<3>>(&surrounding_points[0],
                                                    surrounding_points.size()),
                                ArrayView<double>(&weights[0], weights.size()));

  for (unsigned int i = 0; i < surrounding_points.size(); ++i)
    {
      deallog << "Surrunding point " << i << ' ' << surrounding_points[i]
              << " weight " << weights[i] << std::endl;
    }
  deallog << "Mean point " << new_point_flat << std::endl;
  deallog << "Projected point " << ' ' << new_point << std::endl;
}
