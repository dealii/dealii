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

// Check that we get the correct type of manifold when we clone it. This
// checks that the renamed manifolds are not just aliases but, when cloned,
// give the correct object type.

#include <deal.II/opencascade/manifold_lib.h>

#include <boost/core/demangle.hpp>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepFill.hxx>
#include <GC_MakeCircle.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

#include <typeinfo>

#include "../tests.h"

int
main()
{
  using namespace OpenCASCADE;

  initlog();

  {
    // The circle passing through the vertices of the unit square
    gp_Dir        z_axis(0., 0., 1.);
    gp_Pnt        center(.5, .5, 0.);
    gp_Ax2        axis(center, z_axis);
    Standard_Real radius(std::sqrt(2.) / 2.);

    GC_MakeCircle make_circle(axis, radius);
    Handle(Geom_Circle) circle = make_circle.Value();
    TopoDS_Edge edge           = BRepBuilderAPI_MakeEdge(circle);

    NormalProjectionManifold<2, 3>  manifold_m(edge);
    std::unique_ptr<Manifold<2, 3>> clone_m       = manifold_m.clone();
    const auto                     &deref_clone_m = *clone_m;
    deallog << "typeid of NormalProjectionManifold<2, 3> is "
            << boost::core::demangle(typeid(deref_clone_m).name()) << std::endl;
  }

  {
    // Create a bspline passing through the points
    std::vector<Point<3>> pts1, pts2;
    pts1.push_back(Point<3>(0, 0, 0));
    pts1.push_back(Point<3>(1, 0, 0));

    pts2.push_back(Point<3>(0, 1, 0));
    pts2.push_back(Point<3>(.5, 1, 1));
    pts2.push_back(Point<3>(1, 1, 0));

    TopoDS_Edge edge1 = interpolation_curve(pts1);
    TopoDS_Edge edge2 = interpolation_curve(pts2);

    TopoDS_Face face = BRepFill::Face(edge1, edge2);

    DirectionalProjectionManifold<2, 3> manifold_m(face, Point<3>(0, 0, 1));
    std::unique_ptr<Manifold<2, 3>>     clone_m       = manifold_m.clone();
    const auto                         &deref_clone_m = *clone_m;
    deallog << "typeid of DirectionalProjectionManifold<2, 3> is "
            << boost::core::demangle(typeid(deref_clone_m).name()) << std::endl;
  }

  // Create a bspline passing through the points
  {
    std::vector<Point<3>> pts;
    pts.push_back(Point<3>(0, 0, 0));
    pts.push_back(Point<3>(0, 1, 0));
    pts.push_back(Point<3>(1, 1, 0));
    pts.push_back(Point<3>(1, 0, 0));

    TopoDS_Edge edge1 = interpolation_curve(pts);
    for (unsigned int i = 0; i < pts.size(); ++i)
      pts[i] += Point<3>(0, 0, 1);
    TopoDS_Edge edge2 = interpolation_curve(pts);

    TopoDS_Face face = BRepFill::Face(edge1, edge2);

    NormalToMeshProjectionManifold<1, 3> manifold_m(face);
    std::unique_ptr<Manifold<1, 3>>      clone_m       = manifold_m.clone();
    const auto                          &deref_clone_m = *clone_m;
    deallog << "typeid of NormalProjectionManifold<2, 3> is "
            << boost::core::demangle(typeid(deref_clone_m).name()) << std::endl;
  }

  deallog << "OK" << std::endl;
}
