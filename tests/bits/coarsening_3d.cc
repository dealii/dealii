//----------------------------  coarsening_3d.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  coarsening_3d.cc  ---------------------------


// this test failed with an internal error somewhere in the coarsening
// functions


#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_reordering.h>
#include <grid/grid_generator.h>

#include <fstream>



void create_coarse_grid (Triangulation<3> &coarse_grid)
{
  std::vector<Point<3> >    vertices;
  std::vector<CellData<3> > cells;
  SubCellData               sub_cell_data;
  
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
  coarse_grid.create_triangulation (vertices, cells, sub_cell_data);
}


int main () 
{
  std::ofstream logfile("coarsening_3d.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  Triangulation<3> coarse_grid;
  create_coarse_grid (coarse_grid);

                                   // refine once, then unrefine again
  coarse_grid.refine_global (1);
  for (Triangulation<3>::active_cell_iterator c=coarse_grid.begin_active();
       c != coarse_grid.end(); ++c)
    c->set_coarsen_flag ();
  coarse_grid.execute_coarsening_and_refinement ();
}

  
  
