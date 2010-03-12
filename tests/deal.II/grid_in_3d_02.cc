//----------------------------  grid_in_3d_02.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_in_3d_02.cc  ---------------------------


// check that we can read in and edge orient meshes in 3d with a
// significant number of cells. this tests Mike's reorienting
// algorithm. while the grid_in_3d test checks artificial meshes, this
// one uses a mesh of the 747 nacelle

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
#include <iomanip>

std::ofstream logfile("grid_in_3d_02/output");


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
      gi.read_ucd (in);
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

  GridOut grid_out;
  std::ofstream gnufile("square.gnuplot");
  grid_out.write_gnuplot (tria, gnufile);
}


int main ()
{
  deallog << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ("grid_in_3d_02/747.ucd");
}

