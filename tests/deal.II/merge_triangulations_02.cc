// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
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


// GridGenerator::merge_triangulation did not call
// GridReordering::reorder_cells, but even then it may sometimes fail.
//
// testcase by Carlos Galeano

#include "../tests.h"
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/logstream.h>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iomanip>
#include <cstdio>

std::ofstream logfile("output");

template<int dim>
void mesh_info(const Triangulation<dim> &tria)
{
  deallog << "Mesh info:" << std::endl
	  << " dimension: " << dim << std::endl
	  << " no. of cells: " << tria.n_active_cells() << std::endl;

  // Next loop over all faces of all cells and find how often each boundary
  // indicator is used:
  {
    std::map<unsigned int, unsigned int> boundary_count;
    typename Triangulation<dim>::active_cell_iterator
      cell = tria.begin_active(),
      endc = tria.end();
    for (; cell!=endc; ++cell)
      {
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
          {
            if (cell->face(face)->at_boundary())
              boundary_count[cell->face(face)->boundary_indicator()]++;
          }
      }

    deallog << " boundary indicators: ";
    for (std::map<unsigned int, unsigned int>::iterator it=boundary_count.begin();
         it!=boundary_count.end();
         ++it)
      {
        deallog << it->first << "(" << it->second << " times) ";
      }
    deallog << std::endl;
  }

  // Finally, produce a graphical representation of the mesh to an output
  // file:
  GridOut grid_out;
  grid_out.write_gnuplot (tria, deallog.get_file_stream());
}


void make_grid ()
{
  Triangulation<2> tria1;
  GridGenerator::hyper_cube_with_cylindrical_hole (tria1, 0.25, 1.0);
  
  Triangulation<2> tria3;
  GridGenerator::hyper_cube_with_cylindrical_hole (tria3, 0.25, 1.0);
  GridTools::shift (Point<2>(0,-2), tria3);
  Triangulation<2> triangulation2;

  mesh_info (tria1);
  mesh_info (tria3);
  GridGenerator::merge_triangulations (tria1, tria3, triangulation2);

  mesh_info(triangulation2);
  deallog << "Number of active cells: "
	  << triangulation2.n_active_cells()
	  << std::endl;
  deallog << "Total number of cells: "
	  << triangulation2.n_cells()
	  << std::endl;

}


int main ()
{
  deallog << std::setprecision(2);
  logfile << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  make_grid();
}
