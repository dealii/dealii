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

// Create a sphere, intersect it with a plane, and use the resulting
// object as a boundary.

#include "../tests.h"

#include <deal.II/opencascade/utilities.h>
#include <deal.II/opencascade/boundary_lib.h>

#include <fstream>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_Ax2.hxx>
#include <GC_MakeCircle.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS.hxx>
using namespace OpenCASCADE;

int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  // Create a sphere
  gp_Pnt center(.5,.5,.5);
  Standard_Real radius(Point<3>().distance(point(center)));

  TopoDS_Face face = BRepPrimAPI_MakeSphere(center, radius);

  // Intersect the sphere with the plane z = 0
  TopoDS_Shape intersect = intersect_plane(face, 0, 0, 1, 0);
  // This shape should be the standard circle. Only one edge should be
  // in intersect: our edge.
  TopExp_Explorer exp(intersect, TopAbs_EDGE);
  TopoDS_Edge edge = TopoDS::Edge(exp.Current());

  ArclengthProjectionLineManifold<2,3> boundary_line(edge);

  // The unit square.
  Triangulation<2,3> tria;
  GridGenerator::hyper_cube(tria);

  // Set the exterior boundary
  tria.set_manifold(1, boundary_line);

  // This is here to ignore the points created in the interior of the
  // face.
  tria.begin()->set_all_manifold_ids(1);
  tria.begin()->set_manifold_id(0);

  // We refine twice, and expect the outer points to end up on a
  // smooth curve interpolating the square vertices.
  tria.refine_global(2);


  // You can open the generated file with gmsh.
  GridOut gridout;
  gridout.write_msh (tria, logfile);
  return 0;
}
