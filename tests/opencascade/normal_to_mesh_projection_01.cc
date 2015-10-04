//-----------------------------------------------------------
//
//    Copyright (C) 2014 - 2015 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

// Create a BSpline surface, and test axis projection.

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>

#include <deal.II/opencascade/utilities.h>
#include <deal.II/opencascade/boundary_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_generator.h>

#include <TopTools.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <BRepFill.hxx>
#include <Standard_Stream.hxx>

using namespace OpenCASCADE;

int main ()
{
  std::ofstream logfile("output");

  // Create a bspline passing through the points
  std::vector<Point<3> > pts;
  pts.push_back(Point<3>(0,0,0));
  pts.push_back(Point<3>(0,1,0));
  pts.push_back(Point<3>(1,1,0));
  pts.push_back(Point<3>(1,0,0));

  TopoDS_Edge edge1 = interpolation_curve(pts);
  for (unsigned int i=0; i<pts.size(); ++i)
    pts[i] += Point<3>(0,0,1);
  TopoDS_Edge edge2 = interpolation_curve(pts);

  TopoDS_Face face = BRepFill::Face (edge1, edge2);


  NormalToMeshProjectionBoundary<1,3> manifold(face);

  Triangulation<1,3> tria;
  GridGenerator::hyper_cube(tria);

  tria.begin_active()->set_all_manifold_ids(1);
  tria.set_manifold(1, manifold);

  tria.refine_global(4);
  GridOut gridout;
  gridout.write_msh(tria, logfile);

  return 0;
}
