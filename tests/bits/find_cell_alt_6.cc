// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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



// On a 2D mesh of the following structure look for the cells surrounding
// each vertex, using the find_active_cell_around_point with Mapping:
//
// x-----x-----x
// |     |     |
// |     |     |
// |     |     |
// x--x--x-----x
// |  |  |     |
// x--x--x     x
// |  |  |     |
// x--x--x-----x

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/fe/mapping_q1.h>

#include <fstream>




void check (Triangulation<2> &tria)
{
  const std::vector<Point<2> > &v = tria.get_vertices();
  MappingQ1<2> map;

  for (unsigned i=0; i<tria.n_vertices(); i++)
    {
      std::pair<Triangulation<2>::active_cell_iterator, Point<2> >
      cell = GridTools::find_active_cell_around_point(map, tria, v[i]);

      deallog << "Vertex <" << v[i] << "> found in cell ";
      for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
        deallog << "<" << cell.first->vertex(v) << "> ";
      deallog << " [local: " << cell.second << "]" << std::endl;
    }
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Triangulation<2> coarse_grid;
      GridGenerator::hyper_cube (coarse_grid);
      coarse_grid.refine_global (1);
      coarse_grid.begin_active()->set_refine_flag();
      coarse_grid.execute_coarsening_and_refinement();
      check (coarse_grid);
    }
  catch (const std::exception &exc)
    {
      // we shouldn't get here...
      deallog << "Caught an error..." << std::endl;
      deallog << exc.what() << std::endl;
    }
}



