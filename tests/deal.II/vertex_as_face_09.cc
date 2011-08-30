//----------------------------  vertex_as_face_09.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  vertex_as_face_09.cc  ---------------------------

// verify that we can do things like cell->face() in 1d as well. here:
// Triangulation::get_boundary_indicators() should return an empty vector when
// called to create a mesh that is a loop


#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <fstream>


void test ()
{
  Triangulation<2,2> volume_mesh;
  GridGenerator::hyper_cube (volume_mesh);

  Triangulation<1,2> tria;

  GridTools::extract_boundary_mesh (volume_mesh, tria);

  deallog << "n_cells = " << tria.n_active_cells() << std::endl;
  deallog << "n_boundary_ids = " << tria.get_boundary_indicators ().size()
	  << std::endl;
}



int main ()
{
  std::ofstream logfile("vertex_as_face_09/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test ();

  return 0;
}
