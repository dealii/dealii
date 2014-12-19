// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_reordering.h>

#include <fstream>
#include <iomanip>


// generate two degenerate quads which reduce to the form of triangles, where
// two sides of the quad form one side of the triangle. the neighboring quads
// are neighboring along both of these sides.
//
//  *-------*
//  |1     /|
//  |    a/ |
//  |    /  |
//  |   *   |
//  |  /    |
//  | /b    |
//  |/     2|
//  *-------*
//
// looking from cell 1 at face a and asking for the corresponding face (number)
// of cell 2 neighbor_of_neighbor returns something corresponding to face b,
// not a, as it should. Simply looking at the identity of the neighboring cell
// is not enough, we have to look at the (index of the) face instead.

void create_grid (Triangulation<2> &tria)
{
  const unsigned int n_points=5;

  const Point<2> points[n_points] = { Point<2>(0.0,0.0),
                                      Point<2>(1.0,1.0),
                                      Point<2>(2.0,2.0),
                                      Point<2>(0.0,2.0),
                                      Point<2>(2.0,0.0)
                                    };
  std::vector<Point<2> > vertices(n_points);
  vertices.assign(points,points+n_points);

  std::vector<CellData<2> > cells(2);
  cells[0].vertices[0]   = 0;
  cells[0].vertices[1]   = 1;
  cells[0].vertices[2]   = 2;
  cells[0].vertices[3]   = 3;
  cells[1].vertices[0]   = 0;
  cells[1].vertices[1]   = 4;
  cells[1].vertices[2]   = 2;
  cells[1].vertices[3]   = 1;

  // generate a triangulation
  // out of this
  GridReordering<2>::reorder_cells (cells);
  try
    {
      tria.create_triangulation_compatibility (vertices, cells, SubCellData());
    }
  catch (Triangulation<2>::DistortedCellList &dcv)
    {
      // ignore the exception
    }
}


void check_neighbors (const Triangulation<2> &tria)
{
  Triangulation<2>::cell_iterator cell = tria.begin();
  for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
    if (cell->neighbor(f).state()==IteratorState::valid)
      {
        const unsigned int neighbor_neighbor=cell->neighbor_of_neighbor(f);
        deallog << "At face " << f << ": neighbor_of_neighbor="
                << neighbor_neighbor << std::endl;
        Assert(cell->face(f)==cell->neighbor(f)->face(neighbor_neighbor),
               ExcMessage("Error in neighbor_of_neighbor() function!"));
      }
}



int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Triangulation<2> tria;
  create_grid(tria);
  check_neighbors(tria);
}

