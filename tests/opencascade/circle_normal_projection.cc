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

// Create a circle, a Triangulation, and try to project normally on
// it.

#include "../tests.h"

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
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

using namespace OpenCASCADE;

int main () 
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

				   // The circle passing through the
				   // vertices of the unit square
  gp_Dir z_axis(0.,0.,1.);
  gp_Pnt center(.5,.5,0.);
  gp_Ax2 axis(center, z_axis);
  Standard_Real radius(std::sqrt(2.)/2.);
  
  Handle(Geom_Curve) circle = GC_MakeCircle(axis, radius);
  TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(circle);
  
  
				   // Create a boundary projector.
  NormalProjectionBoundary<2,3> boundary_line(edge);
  
				   // This one is for checking: This
				   // is what deal.II would do for a
				   // circle.
  HyperBallBoundary<2,3> boundary_line_deal (Point<3>(.5,.5,0),
					     std::sqrt(2.)/2.);
  
  
  
				   // The unit square.
  Triangulation<2,3> tria;
  GridGenerator::hyper_cube(tria);

				   // Set the exterior boundary
  tria.set_boundary(0, boundary_line);

				   // This is here to ignore the
				   // points created in the interior
				   // of the face.
  tria.begin()->set_material_id(1);

				   // We refine twice, and expect the
				   // outer points to end up on the
				   // circle.
  tria.refine_global(2);


				   // You can open the generated file
				   // with paraview.
  GridOut gridout;
  gridout.write_ucd (tria, logfile);
  
  return 0;
}
                  
