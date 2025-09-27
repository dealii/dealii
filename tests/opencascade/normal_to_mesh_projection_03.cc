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


// Create a BSpline surface, and test axis projection for Q2 elements.

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1.h>

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
  pts.push_back(Point<3>(1, 0, -1));
  pts.push_back(Point<3>(0, 0, -1));

  TopoDS_Edge edge1 = interpolation_curve(pts);
  for (unsigned int i = 0; i < pts.size(); ++i)
    pts[i] += Point<3>(0, 1, 0);
  TopoDS_Edge edge2 = interpolation_curve(pts);

  TopoDS_Face face = BRepFill::Face(edge1, edge2);


  NormalToMeshProjectionManifold<2, 3> manifold(face);

  Triangulation<2, 3> tria;
  GridGenerator::hyper_cube(tria);

  tria.begin_active()->set_all_manifold_ids(1);
  tria.set_manifold(1, manifold);

  tria.refine_global(1);

  Tensor<1, 3> direction({.1, .01, .001});
  double       tolerance = 1e-10;

  {
    FE_Q<2, 3>       fe(2);
    DoFHandler<2, 3> dh(tria);
    dh.distribute_dofs(fe);
    MappingQ<2, 3>        mapping2(2);
    std::vector<Point<3>> spoints(dh.n_dofs());
    DoFTools::map_dofs_to_support_points(mapping2, dh, spoints);

    std::sort(spoints.begin(),
              spoints.end(),
              [&](const Point<3> &p1, const Point<3> &p2) {
                return OpenCASCADE::point_compare(p1, p2, direction, tolerance);
              });

    for (auto p : spoints)
      deallog << p << std::endl;
    deallog << "============" << std::endl;
  }

  tria.refine_global(1);
  {
    FE_Q<2, 3>       fe(1);
    DoFHandler<2, 3> dh(tria);
    dh.distribute_dofs(fe);
    std::vector<Point<3>> spoints(dh.n_dofs());
    DoFTools::map_dofs_to_support_points(StaticMappingQ1<2, 3>::mapping,
                                         dh,
                                         spoints);

    std::sort(spoints.begin(),
              spoints.end(),
              [&](const Point<3> &p1, const Point<3> &p2) {
                return OpenCASCADE::point_compare(p1, p2, direction, tolerance);
              });

    for (auto p : spoints)
      deallog << p << std::endl;
  }
}
