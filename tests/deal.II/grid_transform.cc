//----------------------------  grid_transform.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_transform.cc  ---------------------------



#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <fe/mapping_q.h>

#include <fstream>


int main ()
{
  const unsigned int dim=2;
  Point<dim> origin;
  HyperShellBoundary<dim> boundary(origin);
  MappingQ<dim> mapping(2);
  Triangulation<dim> tria;
  tria.set_boundary(0, boundary);
  const double inner_radius=1.;
  const double outer_radius=5.;
  GridGenerator::hyper_shell(tria, origin, inner_radius, outer_radius, 8);
  tria.refine_global(2);
  
  GridOut grid_out;
  if (false)
    {
      std::ofstream eps_stream1("grid_untransformed.output");
      grid_out.write_eps(tria, eps_stream1, &mapping);
    }
  
				   // build up a map of vertex indices
				   // of boundary vertices to the new
				   // boundary points
  std::map<unsigned int,Point<dim> > new_points;
  
				   // new center and new radius
				   // of the inner circle.
  const Point<dim> n_center(0,-1);
  const double n_radius=0.5;

  Triangulation<dim>::cell_iterator cell=tria.begin_active(),
					     endc=tria.end();
  Triangulation<dim>::face_iterator face;
  for (; cell!=endc; ++cell)
    {
      if (cell->at_boundary())
	for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	  {
	    face=cell->face(face_no);
	    if (face->at_boundary())
	      for (unsigned int vertex_no=0;
		   vertex_no<GeometryInfo<dim>::vertices_per_face; ++vertex_no)
		{
		  const Point<dim> &v=face->vertex(vertex_no);
		  if (fabs(sqrt(v.square())-outer_radius)<1e-15)
		    {
						       // leave the
						       // point, where
						       // they are.
		      new_points.insert(std::pair<unsigned int, Point<dim> > (
			face->vertex_index(vertex_no), v));
		    }
		  else if (fabs(sqrt(v.square())-inner_radius)<1e-15)
		    {
						       // move the
						       // center of
						       // the inner
						       // circle to
						       // (-1,0) and
						       // take half
						       // the radius
						       // of the
						       // circle.
		      new_points.insert(std::pair<unsigned int, Point<dim> > (
			face->vertex_index(vertex_no), n_radius/inner_radius*v+n_center));
		      face->set_boundary_indicator(1);
		    }
		  else
		    Assert(false, ExcInternalError());
		}  
	  }
    }

				   // test output
  if (false)
    {
      std::cout << "New points:" << std::endl;
      std::map<unsigned int,Point<dim> >::const_iterator map_iter;
      for (map_iter=new_points.begin(); map_iter!=new_points.end(); ++map_iter)
	std::cout << "vertex " << map_iter->first << ": " << map_iter->second << std::endl;
    }

  GridGenerator::laplace_transformation (tria, new_points);
  HyperBallBoundary<dim> inner_ball(n_center, n_radius);
  tria.set_boundary(1, inner_ball);
  std::ofstream eps_stream2("grid_transform.output");
  grid_out.write_eps(tria, eps_stream2, &mapping);

  tria.clear();
  
  return 0;
};
