//----------------------------  geometry_info_2.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  geometry_info_2.cc  ---------------------------


// output all integer values and functions of GeometryInfo

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/geometry_info.h>

#include <fstream>
#include <cstdlib>

template <int dim>
void test ()
{
  deallog << "max_children_per_cell "
	  << GeometryInfo<dim>::max_children_per_cell << std::endl;
  deallog << "faces_per_cell    "
	  << GeometryInfo<dim>::faces_per_cell << std::endl;
  deallog << "max_children_per_face "
	  << GeometryInfo<dim>::max_children_per_face << std::endl;
  deallog << "vertices_per_cell "
	  <<GeometryInfo<dim>::vertices_per_cell << std::endl;
  deallog << "lines_per_cell    "
	  << GeometryInfo<dim>::lines_per_cell << std::endl;
  deallog << "quads_per_cell    "
	  << GeometryInfo<dim>::quads_per_cell << std::endl;
  deallog << "hexes_per_cell    "
	  << GeometryInfo<dim>::hexes_per_cell << std::endl;

  deallog << "vertices_per_face "
	  << GeometryInfo<dim>::vertices_per_face << std::endl;
  deallog << "lines_per_face    "
	  << GeometryInfo<dim>::lines_per_face << std::endl;
  deallog << "quads_per_face    "
	  << GeometryInfo<dim>::quads_per_face << std::endl;

  for (unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell;++f)
    deallog << "face normal" << f << ' '
	    << (GeometryInfo<dim>::unit_normal_orientation[f] > 0.
		? '+' : '-')
	    << "x" << GeometryInfo<dim>::unit_normal_direction[f]
	    << std::endl;

  for (unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell;++f)
    {
      deallog << "face_children" << f << "[true ]";
      for (unsigned int v=0;v < GeometryInfo<dim>::max_children_per_face;++v)
	deallog << ' '
		<< GeometryInfo<dim>::child_cell_on_face(RefinementCase<dim>::isotropic_refinement,f, v, true);
      deallog << std::endl;
      deallog << "face_children" << f << "[false]";
      for (unsigned int v=0;v < GeometryInfo<dim>::max_children_per_face;++v)
	deallog << ' '
		<< GeometryInfo<dim>::child_cell_on_face(RefinementCase<dim>::isotropic_refinement,f, v, false);
      deallog << std::endl;
    }
  
  for (unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell;++f)
    {
      deallog << "face_vertices" << f << "[true ]";
      for (unsigned int v=0;v < GeometryInfo<dim>::vertices_per_face;++v)
	deallog << ' '
		<< GeometryInfo<dim>::face_to_cell_vertices(f, v, true);
      deallog << std::endl;
      deallog << "face_vertices" << f << "[false]";
      for (unsigned int v=0;v < GeometryInfo<dim>::vertices_per_face;++v)
	deallog << ' '
		<< GeometryInfo<dim>::face_to_cell_vertices(f, v, false);
      deallog << std::endl;
    }
  
  for (unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell;++f)
    {
      deallog << "face_lines" << f << "[true ]";
      for (unsigned int v=1;v <= GeometryInfo<dim>::lines_per_face;++v)
	deallog << ' '
		<< GeometryInfo<dim>::face_to_cell_lines(f, v-1, true);
      deallog << std::endl;
      deallog << "face_lines" << f << "[false]";
      for (unsigned int v=1;v <= GeometryInfo<dim>::lines_per_face;++v)
	deallog << ' '
		<< GeometryInfo<dim>::face_to_cell_lines(f, v-1, false);
      deallog << std::endl;
    }
  
  for (unsigned int f=0;f<GeometryInfo<dim>::lines_per_cell;++f)
    {
      deallog << "line_vertices" << f;
      for (unsigned int v=0;v < GeometryInfo<1>::vertices_per_cell;++v)
	deallog << ' '
		<< GeometryInfo<dim>::line_to_cell_vertices(f, v);
      deallog << std::endl;
    }  
}


int main () 
{
  std::ofstream logfile("geometry_info_2/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("1d");
  test<1> ();
  deallog.pop();
  deallog.push("2d");
  test<2> ();
  deallog.pop();
  deallog.push("3d");
  test<3> ();
  deallog.pop();
  return 0;
}
