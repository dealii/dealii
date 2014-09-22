// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// check that we can read in and edge orient meshes in 3d with a
// significant number of cells. this tests Mike's reorienting
// algorithm

#include "../tests.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <iomanip>

std::ofstream logfile("output");


void test (const char *filename)
{
  const unsigned int dim = 3;

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

void test1()
{
  Triangulation<3> tria;
  GridIn<3> gi;
  gi.attach_triangulation (tria);
  gi.read ("1.nc");

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

  test (SOURCE_DIR "/grid_in_3d/1.in");
  test (SOURCE_DIR "/grid_in_3d/2.in");
  test (SOURCE_DIR "/grid_in_3d/3.in");
  test (SOURCE_DIR "/grid_in_3d/4.in");

  test (SOURCE_DIR "/grid_in_3d/evil_0.in");
  test (SOURCE_DIR "/grid_in_3d/evil_1.in");
  test (SOURCE_DIR "/grid_in_3d/evil_2.in");
  test (SOURCE_DIR "/grid_in_3d/evil_3.in");
  test (SOURCE_DIR "/grid_in_3d/evil_4.in");

  // test1 needs NetCDF
//    test1 ();
}

