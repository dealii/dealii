// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// A simplified version of the _08 test that uses fewer cells

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



void
cylinder(Triangulation<3> &tria,
         const double      radius      = 1,
         const double      half_length = 1)
{
  // Copy the base from hyper_ball<3>
  // and transform it to yz
  const double d            = radius / std::sqrt(2.0);
  const double a            = d / (1 + std::sqrt(2.0));
  Point<3>     vertices[16] = {
    Point<3>(-d, -half_length, -d),
    Point<3>(d, -half_length, -d),
    Point<3>(-a, -half_length, -a),
    Point<3>(a, -half_length, -a),
    Point<3>(-a, -half_length, a),
    Point<3>(a, -half_length, a),
    Point<3>(-d, -half_length, d),
    Point<3>(d, -half_length, d),
    Point<3>(-d, 0, -d),
    Point<3>(d, 0, -d),
    Point<3>(-a, 0, -a),
    Point<3>(a, 0, -a),
    Point<3>(-a, 0, a),
    Point<3>(a, 0, a),
    Point<3>(-d, 0, d),
    Point<3>(d, 0, d),
  };
  // Turn cylinder such that y->x
  for (unsigned int i = 0; i < 16; ++i)
    {
      const double h = vertices[i][1];
      vertices[i][1] = -vertices[i][0];
      vertices[i][0] = h;
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

  tria.create_triangulation(std::vector<Point<3>>(&vertices[0], &vertices[16]),
                            cells,
                            SubCellData()); // no boundary information
}


void
test()
{
  const int dim = 3;

  Triangulation<dim> triangulation;
  cylinder(triangulation);

  GridOut().write_gnuplot(triangulation, deallog.get_file_stream());

  Triangulation<dim - 1, dim> triangulation_surface;
  GridGenerator::extract_boundary_mesh(triangulation, triangulation_surface);

  GridOut().write_gnuplot(triangulation_surface, deallog.get_file_stream());

  deallog << triangulation_surface.n_used_vertices() << std::endl;
  deallog << triangulation_surface.n_active_cells() << std::endl;
}


int
main()
{
  initlog();

  test();
}
