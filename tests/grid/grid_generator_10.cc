// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// Test GridReordering<3>::reorder_cells

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



void my_cylinder(Triangulation<3> &tria,
                 const double      radius,
                 const double      half_length)
{
  // Copy the base from hyper_ball<3>
  // and transform it to yz
  const double d            = radius / std::sqrt(2.0);
  const double a            = d / (1 + std::sqrt(2.0));
  Point<3>     vertices[16] = {
    Point<3>(-d, 0, -d),
    Point<3>(d, 0, -d),
    Point<3>(-a, 0, -a),
    Point<3>(a, 0, -a),
    Point<3>(-a, 0, a),
    Point<3>(a, 0, a),
    Point<3>(-d, 0, d),
    Point<3>(d, 0, d),
    Point<3>(-d, half_length, -d),
    Point<3>(d, half_length, -d),
    Point<3>(-a, half_length, -a),
    Point<3>(a, half_length, -a),
    Point<3>(-a, half_length, a),
    Point<3>(a, half_length, a),
    Point<3>(-d, half_length, d),
    Point<3>(d, half_length, d),
  };
  // Turn cylinder such that y->x
  for (unsigned int i = 0; i < 16; ++i)
    {
      const double h = vertices[i](1);
      vertices[i](1) = -vertices[i](0);
      vertices[i](0) = h;
    }

  int cell_vertices[5][8] = {{0, 1, 8, 9, 2, 3, 10, 11},
                             {0, 2, 8, 10, 6, 4, 14, 12},
                             {2, 3, 10, 11, 4, 5, 12, 13},
                             {1, 7, 9, 15, 3, 5, 11, 13},
                             {6, 4, 14, 12, 7, 5, 15, 13}};

  std::vector<CellData<3>> cells(5, CellData<3>());

  for (unsigned int i = 0; i < 5; ++i)
    {
      for (unsigned int j = 0; j < 8; ++j)
        cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };

  GridReordering<3>::reorder_cells(cells);
  tria.create_triangulation(std::vector<Point<3>>(&vertices[0], &vertices[16]),
                            cells,
                            SubCellData()); // no boundary information
}


void
check_grid()
{
  const unsigned int dim = 3;
  Triangulation<dim> triangulation;

  my_cylinder(triangulation, 0.5, 1.0);

  Triangulation<dim>::active_cell_iterator cell = triangulation.begin(),
                                           endc = triangulation.end();
  for (; cell != endc; ++cell)
    {
      deallog << cell << std::endl;
      for (const unsigned int face : GeometryInfo<dim>::face_indices())
        {
          deallog << face << ": "
                  << (cell->face_orientation(face) ? "true " : "false ")
                  << (cell->face_flip(face) ? "true " : "false ")
                  << (cell->face_rotation(face) ? "true" : "false")
                  << std::endl;
        }
      for (unsigned int line = 0; line < GeometryInfo<dim>::lines_per_cell;
           ++line)
        {
          deallog << line << ": "
                  << (cell->line_orientation(line) ? "true" : "false")
                  << std::endl;
          Assert(cell->line_orientation(line) == true, ExcInternalError());
        }
    }
}


int
main()
{
  initlog();
  deallog.depth_console(0);

  check_grid();
}
