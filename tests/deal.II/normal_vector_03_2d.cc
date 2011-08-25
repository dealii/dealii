//----------------------------  normal_vector_03_2d.cc  ---------------------------
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
//----------------------------  normal_vector_03_2d.cc  ---------------------------


// like _01 but instead for random points on the face. note that the normal
// vector on a bilinearly mapped face is also a bilinear function



#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iomanip>



void create_triangulation(const unsigned int case_no,
			  Triangulation<2> &tria)
{
  switch (case_no)
    {
      case 0:
	    GridGenerator::hyper_cube(tria, 1., 3.);
	    break;
      case 1:
      {
	GridGenerator::hyper_cube(tria, 1., 3.);
	Point<2> &v0=tria.begin_active()->vertex(0);
	v0 = Point<2> (-0.5,-1);
	Point<2> &v1=tria.begin_active()->vertex(1);
	v1 = Point<2> (1.25, 0.25);
	break;
      }
      default:
	    Assert(false, ExcNotImplemented());
    };
}



int main ()
{
  std::ofstream logfile ("normal_vector_03_2d/output");
  deallog << std::setprecision (3);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  FE_Q<2> linear_interpolator(1);

  Triangulation<2> tria;
  StraightBoundary<2> boundary;
  Boundary<2>::FaceVertexNormals normals;
  for (unsigned int case_no=0; case_no<2; ++case_no)
    {
      deallog << "Case" << case_no << std::endl;
      create_triangulation(case_no, tria);
      const Triangulation<2>::active_cell_iterator cell=tria.begin_active();
      Triangulation<2>::face_iterator face;
      for (unsigned int face_no=0; face_no<GeometryInfo<2>::faces_per_cell; ++face_no)
	{
	  face=cell->face(face_no);
	  boundary.get_normals_at_vertices(face, normals);

	  for (double xi=0; xi<=1; xi+=0.234)
	    for (double eta=0; eta<=1; eta+=0.234)
	      {
		Point<2> p;
		Tensor<1,2> normal;

		for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_face; ++v)
		  {
		    p += face->vertex(v) * linear_interpolator.shape_value(v,Point<2>(xi,eta));
		    normal += normals[v] *
			      linear_interpolator.shape_value(v,Point<2>(xi,eta));
		  }
		normal /= normal.norm();

		deallog << "p=" << p
			<< ", n=" << boundary.normal_vector (face, p)
			<< std::endl;

		Assert ((boundary.normal_vector (face, p)
			 -
			normal).norm()
			<
			1e-10,
			ExcInternalError());
	      }
	}

      tria.clear();
    }

  deallog << "OK" << std::endl;
}
