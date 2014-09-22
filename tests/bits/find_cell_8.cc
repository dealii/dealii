// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// take a 2d mesh and check that we can find an arbitrary point's cell
// in it. We consider a special triangulation, where the point p does not lie
// in a cell adjacent to the vertex with minimal distance to p. The test should
// fail for all revisions <= 25704M.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <fstream>

void create_coarse_grid(Triangulation<2> &coarse_grid)
{

  static const Point<2> vertices_1[]
    = {  Point<2> (0.,  0.),//0
         Point<2> (1.,  0.),//1
         Point<2> (1.,  1.),//2
         Point<2> (0.,  1.),//3
         Point<2> (1.3,   0.),//4
         Point<2> (1.3, 1./2.),//5
      };
  const unsigned int
  n_vertices = sizeof(vertices_1) / sizeof(vertices_1[0]);

  const std::vector<Point<2> > vertices (&vertices_1[0],
                                         &vertices_1[n_vertices]);

  static const int cell_vertices[][GeometryInfo<2>::vertices_per_cell]
  = {{0, 1, 3, 2},
    {1, 4, 2, 5},
  };
  const unsigned int
  n_cells = sizeof(cell_vertices) / sizeof(cell_vertices[0]);

  std::vector<CellData<2> > cells (n_cells, CellData<2>());
  for (unsigned int i=0; i<n_cells; ++i)
    {
      for (unsigned int j=0;
           j<GeometryInfo<2>::vertices_per_cell;
           ++j)
        cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    }

  coarse_grid.create_triangulation (vertices,
                                    cells,
                                    SubCellData());
}


void check (Triangulation<2> &tria)
{
  Point<2> p (0.99, 1./2.);

  Triangulation<2>::active_cell_iterator cell
    = GridTools::find_active_cell_around_point (tria, p);

  deallog << cell << std::endl;
  for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
    deallog << "<" << cell->vertex(v) << "> ";
  deallog << std::endl;

  Assert (p.distance (cell->center()) < cell->diameter()/2,
          ExcInternalError());
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  {
    Triangulation<2> coarse_grid;
    create_coarse_grid(coarse_grid);
    check (coarse_grid);
  }
}



