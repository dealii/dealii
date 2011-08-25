//----------------------------  normals_at_vertices_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2008, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  normals_at_vertices_01.cc  ---------------------------


/* Author: Ralf Hartmann, 2005 */



#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>

#include <fstream>
#include <iomanip>



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
  std::ofstream logfile ("normals_at_vertices_01/output");
  deallog << std::setprecision (3);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  Triangulation<3> tria;
  StraightBoundary<3> boundary;
  Boundary<3>::FaceVertexNormals normals;
  for (unsigned int case_no=0; case_no<2; ++case_no)
    {
      deallog << "Case" << case_no << std::endl;
      create_triangulation(case_no, tria);
      const Triangulation<3>::active_cell_iterator cell=tria.begin_active();
      Triangulation<3>::face_iterator face;
      for (unsigned int face_no=0; face_no<GeometryInfo<3>::faces_per_cell; ++face_no)
	{
	  deallog << " Face" << face_no << std::endl;
	  face=cell->face(face_no);
	  boundary.get_normals_at_vertices(face, normals);
	  for (unsigned int v=0; v<GeometryInfo<3>::vertices_per_face; ++v)
	    {
	      deallog << "  vertex=" << face->vertex(v)
		      << ",  normal=" << normals[v] << std::endl;

					       // note that we can't check
					       // here that the normal vector
					       // is, in fact, normalized,
					       // since the function does not
					       // actually guarantee that
	    }
	}
      tria.clear();
    }
}
