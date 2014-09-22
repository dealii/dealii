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


// read a file in the MSH format used by the GMSH program. test the
// reader for version 2 of the MSH file format

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
#include <string>

std::ofstream logfile("output");


template<int dim>
void check_file (const std::string name,
                 typename GridIn<dim>::Format format)
{
  Triangulation<dim> tria;
  GridIn<dim> gi;
  gi.attach_triangulation (tria);
  gi.read(name, format);
  deallog << '\t' << tria.n_vertices()
          << '\t' << tria.n_cells()
          << std::endl;

  GridOut grid_out;
  grid_out.write_gnuplot (tria, deallog.get_file_stream());
}

void filename_resolution()
{
  check_file<2> (std::string(SOURCE_DIR "/grid_in_msh_version_2/hole81"), GridIn<2>::msh);
  check_file<2> (std::string(SOURCE_DIR "/grid_in_msh_version_2/hole8170"), GridIn<2>::msh);
}


int main ()
{
  deallog << std::setprecision (5);
  logfile << std::setprecision (5);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  filename_resolution();
}

