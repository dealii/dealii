// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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


// read a 3d file in the VTK format

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

  // now also output boundary ids for all boundary faces. this is
  // where the bug was for which this test was written
  deallog << "Face boundary_ids:" << std::endl;
  for (typename Triangulation<dim>::active_cell_iterator
       cell = tria.begin_active();
       cell != tria.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->at_boundary(f))
        deallog << cell << ' ' << f << ": "
                << (int)cell->face(f)->boundary_id()
                << std::endl;
}

void filename_resolution()
{
  check_file<3> (std::string(SOURCE_DIR "/grid_in_vtk_3d_04.vtk"), GridIn<3>::vtk);
}


int main ()
{
  deallog << std::setprecision (2);
  logfile << std::setprecision (2);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  filename_resolution();
}

