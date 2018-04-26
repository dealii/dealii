//-----------------------------------------------------------
//
//    Copyright (C) 2014 - 2018 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

// Create a BSpline surface, and test axis projection for Q2 elements.

#include "../tests.h"

#include <deal.II/opencascade/utilities.h>
#include <deal.II/opencascade/boundary_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <TopTools.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <BRepFill.hxx>
#include <Standard_Stream.hxx>

using namespace OpenCASCADE;

int main ()
{
  initlog();

  // Create a bspline passing through the points
  std::vector<Point<3> > pts;
  pts.push_back(Point<3>(0,0,0));
  pts.push_back(Point<3>(1,0,0));
  pts.push_back(Point<3>(1,0,-1));
  pts.push_back(Point<3>(0,0,-1));

  TopoDS_Edge edge1 = interpolation_curve(pts);
  for (unsigned int i=0; i<pts.size(); ++i)
    pts[i] += Point<3>(0,1,0);
  TopoDS_Edge edge2 = interpolation_curve(pts);

  TopoDS_Face face = BRepFill::Face (edge1, edge2);


  NormalToMeshProjectionBoundary<2,3> manifold(face);

  Triangulation<2,3> tria;
  GridGenerator::hyper_cube(tria);

  tria.begin_active()->set_all_manifold_ids(1);
  tria.set_manifold(1, manifold);

  tria.refine_global(1);

  Tensor<1,3> direction({.1, .01, .001});
  double tolerance = 1e-10;

  {
    FE_Q<2,3> fe(2);
    DoFHandler<2,3> dh(tria);
    dh.distribute_dofs(fe);
    MappingQGeneric<2,3> mapping2(2);
    std::vector<Point<3>> spoints(dh.n_dofs());
    DoFTools::map_dofs_to_support_points(mapping2, dh, spoints);

    std::sort(spoints.begin(), spoints.end(),
              [&](const Point<3> &p1, const Point<3> &p2)
    {
      return OpenCASCADE::point_compare(p1, p2, direction, tolerance);
    });

    for (auto p: spoints)
      deallog << p << std::endl;
    deallog << "============" << std::endl;
  }

  tria.refine_global(1);
  {
    FE_Q<2,3> fe(1);
    DoFHandler<2,3> dh(tria);
    dh.distribute_dofs(fe);
    std::vector<Point<3>> spoints(dh.n_dofs());
    DoFTools::map_dofs_to_support_points(StaticMappingQ1<2,3>::mapping,
                                         dh, spoints);

    std::sort(spoints.begin(), spoints.end(),
              [&](const Point<3> &p1, const Point<3> &p2)
    {
      return OpenCASCADE::point_compare(p1, p2, direction, tolerance);
    });

    for (auto p: spoints)
      deallog << p << std::endl;
  }
}
