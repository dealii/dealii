//----------------------------  mesh_3d.h  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_3d.h  ---------------------------

// generate two cubes that are attached to each other in a way so that
// the edges are all ok, but the normals of the common face don't
// match up for the standard orientation of the normals. we thus have
// to store the face orientation in each cell

#include "../tests.h"
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_reordering.h>



void create_two_cubes (Triangulation<3> &coarse_grid)
{
  const Point<3> points[6] = { Point<3>(0,0,0),
                               Point<3>(1,0,0),
                               Point<3>(1,1,0),
                               Point<3>(0,1,0),
                               Point<3>(2,0,0),
                               Point<3>(2,1,0)};
  std::vector<Point<3> > vertices;
  for (unsigned int i=0; i<6; ++i)
    vertices.push_back (points[i]);
  for (unsigned int i=0; i<6; ++i)
    vertices.push_back (points[i] + Point<3>(0,0,-1));

  std::vector<CellData<3> > cells;
    
  const unsigned int connectivity[2][4]
    = { { 0,1,2,3 },
        { 4,5,2,1 } };
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
  coarse_grid.create_triangulation (vertices, cells, SubCellData());
}
  




void create_L_shape (Triangulation<3> &coarse_grid)
{
  std::vector<Point<3> >    vertices;
  std::vector<CellData<3> > cells;
  
  const Point<3> outer_points[8] = { Point<3>(-1,0,0),
                                     Point<3>(-1,-1,0),
                                     Point<3>(0,-1,0),
                                     Point<3>(+1,-1,0),
                                     Point<3>(+1,0,0),
                                     Point<3>(+1,+1,0),
                                     Point<3>(0,+1,0),
                                     Point<3>(-1,+1,0) };

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
        { 0, 5, 6, 7 } };
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
  coarse_grid.create_triangulation (vertices, cells, SubCellData());
}


void coarsen_global (Triangulation<3> &grid)
{
  for (Triangulation<3>::active_cell_iterator c=grid.begin_active();
       c != grid.end(); ++c)
    c->set_coarsen_flag ();
  grid.execute_coarsening_and_refinement ();
}
