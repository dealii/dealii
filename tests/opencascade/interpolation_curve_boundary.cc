
//----------------------------  occ_circle_projection.cc  ---------------------------
//    $Id: testsuite.html 13373 2006-07-13 13:12:08Z kanschat $
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  occ_circle_projection.cc  ---------------------------


// Create a Triangulation, interpolate its boundary points to a close
// smooth BSpline, and use that as a Boundary Descriptor.

#include "../tests.h"

#include <deal.II/grid/occ_boundary_lib.h>

#include <fstream>
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/grid_out.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>

#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_Ax2.hxx>
#include <GC_MakeCircle.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

using namespace OpenCASCADE;

int main () 
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  // The curve passing through the vertices of the unit square.
  std::vector<Point<3> > vertices;
  vertices.push_back(Point<3>(0,0,0));
  vertices.push_back(Point<3>(1,0,0));
  vertices.push_back(Point<3>(1,1,0));
  vertices.push_back(Point<3>(0,1,0));
  TopoDS_Shape shape = interpolation_curve(vertices, Point<3>(), true);
  
  // Create a boundary projector.
  NormalProjectionBoundary<2,3> boundary_line(shape);

  // The unit square.
  Triangulation<2,3> tria;
  GridGenerator::hyper_cube(tria);

  // Set the exterior boundary
  tria.set_boundary(0, boundary_line);

  // This is here to ignore the points created in the interior of the
  // face.
  tria.begin()->set_material_id(1);

  // We refine twice, and expect the outer points to end up on a
  // smooth curve interpolating the square vertices.
  tria.refine_global(2);


  // You can open the generated file with paraview.
  GridOut gridout;
  gridout.write_msh (tria, logfile);

  std::ofstream of("/tmp/grid.msh");
  gridout.write_msh(tria, of);
  
  return 0;
}
                  
