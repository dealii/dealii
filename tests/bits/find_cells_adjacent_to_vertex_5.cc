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



// similar to the _3 test: a mesh where 8 2d cells meet at one
// point. we then extrude this mesh into the z-direction
//
// this test does not use any local refinement

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>

#include <fstream>


void check (Triangulation<3> &tria)
{
  for (unsigned i=0; i<tria.n_vertices(); i++)
    {
      std::vector<Triangulation<3>::active_cell_iterator>
      cells = GridTools::find_cells_adjacent_to_vertex(tria, i);

      deallog << "Vertex " << i << " at "
              << tria.get_vertices()[i] << ": "
              << cells.size() << " cells" << std::endl;

      for (unsigned c=0; c<cells.size(); c++)
        deallog << "   " << cells[c] << std::endl;
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
      // set up the vertices of the
      // cells: the center plus 16
      // points on the perimeter of a
      // circle
      const unsigned int dim = 3;
      std::vector<Point<dim> > vertices;
      vertices.push_back (Point<dim>(0,0,0));

      for (unsigned int i=0; i<16; ++i)
        vertices.push_back (Point<dim>(std::cos(i*2*numbers::PI/16),
                                       std::sin(i*2*numbers::PI/16),
                                       0));

      // then extrude the mesh into
      // z-direction
      for (unsigned int i=0; i<17; ++i)
        vertices.push_back (Point<dim>(vertices[i][0],
                                       vertices[i][1],
                                       1));

      // now create the 8 cells
      std::vector<CellData<dim> > cells;
      for (unsigned int c=0; c<8; ++c)
        {
          CellData<dim> d;
          d.vertices[0] = 0;
          d.vertices[1] = 1+2*c;
          d.vertices[2] = 1+((2*c+2) % 16);
          d.vertices[3] = 1+2*c+1;

          for (unsigned int v=4; v<8; ++v)
            d.vertices[v] = d.vertices[v-4]+17;

          cells.push_back(d);
        }

      Triangulation<dim> coarse_grid;
      coarse_grid.create_triangulation (vertices, cells,
                                        SubCellData());
      check (coarse_grid);
    }
  catch (const std::exception &exc)
    {
      // we shouldn't get here...
      deallog << "Caught an error..." << std::endl;
      deallog << exc.what() << std::endl;
    }
}
