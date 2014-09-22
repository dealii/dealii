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

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_reordering.h>


// generate two cubes that are attached to each other in a way so that
// the edges are all ok, but the normals of the common face don't
// match up for the standard orientation of the normals. we thus have
// to store the face orientation in each cell

void create_two_cubes (Triangulation<3> &coarse_grid)
{
  const Point<3> points[6] = { Point<3>(0,0,0),
                               Point<3>(1,0,0),
                               Point<3>(1,1,0),
                               Point<3>(0,1,0),
                               Point<3>(2,0,0),
                               Point<3>(2,1,0)
                             };
  std::vector<Point<3> > vertices;
  for (unsigned int i=0; i<6; ++i)
    vertices.push_back (points[i]);
  for (unsigned int i=0; i<6; ++i)
    vertices.push_back (points[i] + Point<3>(0,0,-1));

  std::vector<CellData<3> > cells;

  const unsigned int connectivity[2][4]
  = { { 0,1,2,3 },
    { 4,5,2,1 }
  };
  for (unsigned int i=0; i<2; ++i)
    {
      CellData<3> cell;
      for (unsigned int j=0; j<4; ++j)
        {
          cell.vertices[j]   = connectivity[i][j];
          cell.vertices[j+4] = connectivity[i][j]+6;
        }
      cells.push_back (cell);
    }

  // finally generate a triangulation
  // out of this
  GridReordering<3>::reorder_cells (cells);
  coarse_grid.create_triangulation_compatibility (vertices, cells, SubCellData());
}



// generate two cubes that are attached to each other in a way so that
// the edges are not all ok and the common face is rotated. we thus have
// to store the face rotation (and face flip) in each cell

void create_two_cubes_rotation (Triangulation<3> &coarse_grid,
                                const unsigned int n_rotations)
{
  Assert(n_rotations<4, ExcNotImplemented());

  const Point<3> points[6] = { Point<3>(0,0,0),
                               Point<3>(1,0,0),
                               Point<3>(1,1,0),
                               Point<3>(0,1,0),
                               Point<3>(2,0,0),
                               Point<3>(2,1,0)
                             };
  std::vector<Point<3> > vertices;
  for (unsigned int i=0; i<6; ++i)
    vertices.push_back (points[i]);
  for (unsigned int i=0; i<6; ++i)
    vertices.push_back (points[i] + Point<3>(0,0,-1));

  std::vector<CellData<3> > cells(2);

  const unsigned int connectivity[5][8]
  // first row: left cube
  // second row: right cube, standard_orientation
  // third row: right cube, rotated once
  // forth row: right cube, rotated twice
  // fifth row: right cube, rotated three times
  = { { 0,1,2,3,6,7,8,9 },
    { 1,4,5,2,7,10,11,8},
    { 2,5,11,8,1,4,10,7},
    { 8,11,10,7,2,5,4,1},
    { 7,10,4,1,8,11,5,2}
  };
  for (unsigned int j=0; j<8; ++j)
    {
      cells[0].vertices[j]   = connectivity[0][j];
      cells[1].vertices[j]   = connectivity[1+n_rotations][j];
    }
  // finally generate a triangulation
  // out of this
  coarse_grid.create_triangulation_compatibility (vertices, cells, SubCellData());
}



void create_L_shape (Triangulation<3> &coarse_grid)
{
  std::vector<Point<3> >    vertices;
  std::vector<CellData<3> > cells;

  const Point<3> outer_points[7] = { Point<3>(-1,0,0),
                                     Point<3>(-1,-1,0),
                                     Point<3>(0,-1,0),
                                     Point<3>(+1,-1,0),
                                     Point<3>(+1,0,0),
                                     Point<3>(+1,+1,0),
                                     Point<3>(0,+1,0)
                                   };

  // first the point in the middle
  // and the rest of those on the
  // upper surface
  vertices.push_back (Point<3>(0,0,0));
  for (unsigned int i=0; i<7; ++i)
    vertices.push_back (outer_points[i]);

  // then points on lower surface
  vertices.push_back (Point<3>(0,0,-1));
  for (unsigned int i=0; i<7; ++i)
    vertices.push_back (outer_points[i]
                        +
                        Point<3>(0,0,-1));

  const unsigned int n_vertices_per_surface = 8;
  Assert (vertices.size() == n_vertices_per_surface*2,
          ExcInternalError());

  const unsigned int connectivity[3][4]
  = { { 1, 2, 3, 0 },
    { 3, 4, 5, 0 },
    { 0, 5, 6, 7 }
  };
  for (unsigned int i=0; i<3; ++i)
    {
      CellData<3> cell;
      for (unsigned int j=0; j<4; ++j)
        {
          cell.vertices[j]   = connectivity[i][j];
          cell.vertices[j+4] = connectivity[i][j]+n_vertices_per_surface;
        }
      cells.push_back (cell);
    }

  // finally generate a triangulation
  // out of this
  GridReordering<3>::reorder_cells (cells);
  coarse_grid.create_triangulation_compatibility (vertices, cells, SubCellData());
}


void coarsen_global (Triangulation<3> &grid)
{
  for (Triangulation<3>::active_cell_iterator c=grid.begin_active();
       c != grid.end(); ++c)
    c->set_coarsen_flag ();
  grid.execute_coarsening_and_refinement ();
}
