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

#include <fstream>

namespace PhantomGeometry 
{
  const double r1 = 5;
  const double r2 = 10;
  const double dz = 2;
  
  const double r0 = r1/(1.+std::sqrt(2.0));
}


void create_coarse_grid (Triangulation<3> &coarse_grid)
{
  std::vector<Point<3> >    vertices;
  std::vector<CellData<3> > cells;
  SubCellData               sub_cell_data;
  
                                   // first build up the cells of the
                                   // cylinder
  {
                                     // the vertices in each plane of
                                     // the cylinder are located on
                                     // three concentric rings of
                                     // radii r0, r1, and r2,
                                     // respectively. first generate
                                     // these three rings
    const Point<3> ring_points[8] = { Point<3>(-1,0,0),
                                      Point<3>(-1,-1,0) / std::sqrt(2.),
                                      Point<3>(0,-1,0),
                                      Point<3>(+1,-1,0) / std::sqrt(2.),
                                      Point<3>(+1,0,0),
                                      Point<3>(+1,+1,0) / std::sqrt(2.),
                                      Point<3>(0,+1,0),
                                      Point<3>(-1,+1,0) / std::sqrt(2.) };

                                     // first the point in the middle
                                     // and the rest of those on the
                                     // upper surface
    vertices.push_back (Point<3>(0,0,0));
    for (unsigned int i=0; i<8; ++i)
      vertices.push_back (ring_points[i]);

                                     // then points on lower surface
    vertices.push_back (Point<3>(0,0,-PhantomGeometry::dz));
    for (unsigned int i=0; i<8; ++i)
      vertices.push_back (ring_points[i]
                          +
                          Point<3>(0,0,-PhantomGeometry::dz));

    const unsigned int n_vertices_per_surface = 9;
    Assert (vertices.size() == n_vertices_per_surface*2,
            ExcInternalError());
    
                                     // next create cells from these
                                     // vertices. only store the
                                     // vertices of the upper surface,
                                     // the lower ones are the same
                                     // +12
    const unsigned int connectivity[4][4]
      = { { 1, 2, 3, 0 },
          { 3, 4, 5, 0 },
          { 0, 5, 6, 7 },
          { 1, 0, 7, 8 } };
      
                                     // now create cells out of this
    for (unsigned int i=0; i<4; ++i)
      {
        CellData<3> cell;
        for (unsigned int j=0; j<4; ++j)
          {
            cell.vertices[j]   = connectivity[i][j];
            cell.vertices[j+4] = connectivity[i][j]+n_vertices_per_surface;
          }
        cell.material_id = 0;
        cells.push_back (cell);
      }   
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

  
  
