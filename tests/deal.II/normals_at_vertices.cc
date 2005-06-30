//----------------------------  normals_at_vertices.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  normals_at_vertices.cc  ---------------------------


/* Author: Ralf Hartmann, 2005 */



#include "../tests.h"
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_boundary.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>

#include <fstream>



void create_triangulation(const unsigned int case_no,
			  Triangulation<3> &tria)
{
  switch (case_no)
    {
      case 0:
	    GridGenerator::hyper_cube(tria, 1., 3.);
	    break;
      case 1:
      {
	GridGenerator::hyper_cube(tria, 1., 3.);
	Point<3> &v0=tria.begin_active()->vertex(0);
	v0 = Point<3> (0,-0.5,-1);
	Point<3> &v1=tria.begin_active()->vertex(1);
	v1 = Point<3> (1.25, 0.25, 0.25);
	break;
      }
      default:
	    Assert(false, ExcNotImplemented());
    };
}



int main ()
{
  std::ofstream logfile ("normals_at_vertices.output");
  logfile.precision (3);
  logfile.setf(std::ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console (0);

  Triangulation<3> tria;
  StraightBoundary<3> boundary;
  Boundary<3>::FaceVertexNormals normals;
  for (unsigned int case_no=0; case_no<1; ++case_no)
    {
      deallog << "Case" << case_no << std::endl;
      create_triangulation(case_no, tria);
      const Triangulation<3>::active_cell_iterator cell=tria.begin_active();
      Triangulation<3>::face_iterator face;
      for (unsigned int face_no=0; face_no<1/*GeometryInfo<3>::faces_per_cell*/; ++face_no)
	{
	  deallog << " Face" << face_no << std::endl;
	  face=cell->face(face_no);
	  boundary.get_normals_at_vertices(face, normals);
	  for (unsigned int v=0; v<GeometryInfo<3>::vertices_per_face; ++v)
	    deallog << "  vertex=" << face->vertex(v)
		    << ",  normal=" << normals[v] << std::endl;
	}
      tria.clear();
    }
}
