// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// similar to the _3 test: a mesh where 8 2d cells meet at one
// point. we then extrude this mesh into the z-direction
//
// this test does not use any local refinement

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



void
check(Triangulation<3> &tria)
{
  for (unsigned i = 0; i < tria.n_vertices(); ++i)
    {
      std::vector<Triangulation<3>::active_cell_iterator> cells =
        GridTools::find_cells_adjacent_to_vertex(tria, i);

      deallog << "Vertex " << i << " at " << tria.get_vertices()[i] << ": "
              << cells.size() << " cells" << std::endl;

      for (unsigned c = 0; c < cells.size(); ++c)
        deallog << "   " << cells[c] << std::endl;
    }
}


int
main()
{
  initlog();

  try
    {
      // set up the vertices of the
      // cells: the center plus 16
      // points on the perimeter of a
      // circle
      const unsigned int      dim = 3;
      std::vector<Point<dim>> vertices;
      vertices.push_back(Point<dim>(0, 0, 0));

      for (unsigned int i = 0; i < 16; ++i)
        vertices.push_back(Point<dim>(std::cos(i * 2 * numbers::PI / 16),
                                      std::sin(i * 2 * numbers::PI / 16),
                                      0));

      // then extrude the mesh into
      // z-direction
      for (unsigned int i = 0; i < 17; ++i)
        vertices.push_back(Point<dim>(vertices[i][0], vertices[i][1], 1));

      // now create the 8 cells
      std::vector<CellData<dim>> cells;
      for (unsigned int c = 0; c < 8; ++c)
        {
          CellData<dim> d;
          d.vertices[0] = 0;
          d.vertices[1] = 1 + 2 * c;
          d.vertices[2] = 1 + ((2 * c + 2) % 16);
          d.vertices[3] = 1 + 2 * c + 1;

          for (unsigned int v = 4; v < 8; ++v)
            d.vertices[v] = d.vertices[v - 4] + 17;

          cells.push_back(d);
        }

      Triangulation<dim> coarse_grid;
      coarse_grid.create_triangulation(vertices, cells, SubCellData());
      check(coarse_grid);
    }
  catch (const std::exception &exc)
    {
      // we shouldn't get here...
      deallog << "Caught an error..." << std::endl;
      deallog << exc.what() << std::endl;
    }
}
