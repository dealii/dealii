//----------------------------  geometry_info_3.cc  ---------------------------
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
//----------------------------  geometry_info_3.cc  ---------------------------


// check GeometryInfo::face_to_cell_vertices

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/geometry_info.h>

#include <fstream>
#include <cstdlib>


template <int dim>
void test ()
{
  deallog << "Checking in " << dim << "d" << std::endl;
  
  for (unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell;++f)
    for (unsigned int v=0;v<GeometryInfo<dim>::vertices_per_face;++v)
      {
	deallog << "Face " << f << ", vertex=" << v << ": ";
	deallog << GeometryInfo<dim>::face_to_cell_vertices(f,v,true)
		<< std::endl;
      }

  if (dim == 3)
    for (unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell;++f)
      for (unsigned int v=0;v<GeometryInfo<dim>::vertices_per_face;++v)
	{
	  deallog << "Face " << f << ", vertex=" << v
		  << " (reverse orientation): ";
	  deallog << GeometryInfo<dim>::face_to_cell_vertices(f,v,false)
		  << std::endl;
	}
}


int main () 
{
  std::ofstream logfile("geometry_info_3/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  return 0;
}
