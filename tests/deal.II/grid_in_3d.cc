//----------------------------  grid_in_3d.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_in_3d.cc  ---------------------------


// check that we can read in and edge orient meshes in 3d with a
// significant number of cells. this tests Mike's reorienting
// algorithm

#include "../tests.h"
#include <dofs/dof_handler.h>
#include <grid/tria.h>
#include <grid/tria_boundary.h>
#include <grid/tria_boundary_lib.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_out.h>
#include <grid/grid_in.h>
#include <grid/grid_generator.h>
#include <base/logstream.h>

#include <fstream>


std::ofstream logfile("grid_in_3d.output");


void test (const char *filename)
{
  const unsigned int dim = 3;

  deallog << "Reading " << filename << std::endl;
  
  Triangulation<3> tria;
  GridIn<3> gi;
  gi.attach_triangulation (tria);
  std::ifstream in (filename);

  try
    {
      gi.read_xda (in);
    }
  catch (std::exception &exc)
    {
      deallog << "  caught exception:" << std::endl
	      << exc.what()
	      << std::endl;
      return;
    }

  deallog << "  " << tria.n_active_cells() << " active cells" << std::endl;

  int hash = 0;
  int index = 0;
  for (Triangulation<dim>::active_cell_iterator c=tria.begin_active();
       c!=tria.end(); ++c, ++index)
    for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
      hash += (index * i * c->vertex_index(i)) % (tria.n_active_cells()+1);
  deallog << "  hash=" << hash << std::endl;
}


int main ()
{
  logfile.precision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);

  test ("grid_in_3d_1.in");
  test ("grid_in_3d_2.in");
  test ("grid_in_3d_3.in");
  test ("grid_in_3d_4.in");

  test ("grid_in_3d_evil_0.in");
  test ("grid_in_3d_evil_1.in");
  test ("grid_in_3d_evil_2.in");
  test ("grid_in_3d_evil_3.in");
  test ("grid_in_3d_evil_4.in");
}

