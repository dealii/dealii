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
    for (unsigned int ring=0; ring<3; ++ring)
      for (unsigned int i=0; i<8; ++i)
        vertices.push_back (ring_points[i] * (ring == 0 ? PhantomGeometry::r0 :
                                              ring == 1 ? PhantomGeometry::r1 :
                                                          PhantomGeometry::r2));

                                     // then points on lower surface
    vertices.push_back (Point<3>(0,0,-PhantomGeometry::dz));
    for (unsigned int ring=0; ring<3; ++ring)
      for (unsigned int i=0; i<8; ++i)
        vertices.push_back (ring_points[i] * (ring == 0 ? PhantomGeometry::r0 :
                                              ring == 1 ? PhantomGeometry::r1 :
                                                          PhantomGeometry::r2)
                            +
                            Point<3>(0,0,-PhantomGeometry::dz));

    const unsigned int n_vertices_per_surface = 25;
    Assert (vertices.size() == n_vertices_per_surface*2,
            ExcInternalError());
    
                                     // next create cells from these
                                     // vertices. only store the
                                     // vertices of the upper surface,
                                     // the lower ones are the same
                                     // +12
    {
      const unsigned int connectivity[20][4]
        = { { 1, 2, 3, 0 },  // four cells in the center
            { 3, 4, 5, 0 },
            { 0, 5, 6, 7 },
            { 1, 0, 7, 8 },            
          
            { 9, 10, 2, 1 },  // eight cells of inner ring
            { 10, 11, 3, 2 },
            { 11, 12, 4, 3 },
            { 4, 12, 13, 5 },
            { 5, 13, 14, 6 },
            { 6, 14, 15, 7 },
            { 8, 7, 15, 16 },
            { 9, 1, 8, 16 },

            { 17, 18, 10, 9 },  // eight cells of outer ring
            { 18, 19, 11, 10 },
            { 19, 20, 12, 11 },
            { 12, 20, 21, 13 },
            { 13, 21, 22, 14 },
            { 14, 22, 23, 15 },
            { 16, 15, 23, 24 },
            { 17, 9, 16, 24 }   };

                                       // now create cells out of this
      for (unsigned int i=0; i<20; ++i)
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
    
                                     // associate edges and faces on
                                     // the outer boundary with
                                     // boundary indicator 1. do this
                                     // the same way as above, just
                                     // this time with faces (edges
                                     // follow from this
                                     // immediately. some edges are
                                     // duplicated since they belong
                                     // to more than one cell, but
                                     // that doesn't harm us here)
    {
      const unsigned int connectivity[8][2]
        = { { 17,18 }, { 18, 19 }, { 19, 20 }, { 20, 21 },
            { 21,22 }, { 22, 23 }, { 23, 24 }, { 24, 17 }};

      for (unsigned int i=0; i<8; ++i)
        {
          const CellData<2> face = 
            { { connectivity[i][0]+n_vertices_per_surface,
                connectivity[i][1]+n_vertices_per_surface,
                connectivity[i][1],
                connectivity[i][0] },
              1 };
          sub_cell_data.boundary_quads.push_back (face);

          const CellData<1> edges[4] = 
            { { { connectivity[i][0],    connectivity[i][1]    }, 1 },
              { { connectivity[i][0]+n_vertices_per_surface,
                  connectivity[i][1]+n_vertices_per_surface }, 1 },
              { { connectivity[i][0]+n_vertices_per_surface,
                  connectivity[i][0]    }, 1 },
              { { connectivity[i][1]+n_vertices_per_surface,
                  connectivity[i][1]    }, 1 } };
          for (unsigned int i=0; i<4; ++i)
            sub_cell_data.boundary_lines.push_back (edges[i]);
        }
    }
  }

                                   // the second part is setting the
                                   // half-sphere on top of this
  {
                                     // add four cubes to the top of
                                     // the inner four cells, as well
                                     // as 8 to their outside
    {
                                       // mirror the first nine vertices
                                       // above the surface, and scale
                                       // them to a certain distance
                                       // outward
      const double rx = PhantomGeometry::r1 / (1+std::sqrt(3.0));
      for (unsigned int i=0; i<9; ++i)
        {
          Point<3> p (vertices[i][0],
                      vertices[i][1],
                      i == 0 ?
                      1
                      :
                      std::max(std::fabs(vertices[i][0]),
                               std::fabs(vertices[i][1])));
          vertices.push_back (p / std::sqrt(p.square()) * rx);
        }
      Assert (vertices.size() == 59, ExcInternalError());

                                       // same with the next ring of
                                       // vertices, except that they
                                       // go to r1
      for (unsigned int i=9; i<17; ++i)
        {
          Point<3> p (vertices[i][0],
                      vertices[i][1],
                      std::max(std::fabs(vertices[i][0]),
                               std::fabs(vertices[i][1])));
          vertices.push_back (p / std::sqrt(p.square()) *
                              PhantomGeometry::r1);
        }
      Assert (vertices.size() == 67, ExcInternalError());
      
                                       // make 12 cells out of this
      const unsigned int connectivity[12][4]
        = { { 1, 2, 3, 0 },  // four cells in the center
            { 3, 4, 5, 0 },
            { 0, 5, 6, 7 },
            { 1, 0, 7, 8 },

            { 9, 10, 2, 1 },  // eight cells of inner ring
            { 10, 11, 3, 2 },
            { 11, 12, 4, 3 },
            { 4, 12, 13, 5 },
            { 5, 13, 14, 6 },
            { 6, 14, 15, 7 },
            { 8, 7, 15, 16 },
            { 9, 1, 8, 16 },
        };

      for (unsigned int i=0; i<12; ++i)
        {
          CellData<3> cell;
          for (unsigned int j=0; j<4; ++j)
            {
              cell.vertices[j]   = connectivity[i][j]+50;
              cell.vertices[j+4] = connectivity[i][j];
            }
          cell.material_id = 0;
          cells.push_back (cell);
        }
    }

                                     // assign boundary indicators to
                                     // the faces and edges of these
                                     // cells
    {
                                       // these are the numbers of the
                                       // vertices on the top surface
                                       // of the cylinder, with one
                                       // "wrap-around":
      const unsigned int vertices[9] = 
        { 9, 10, 11, 12, 13, 14, 15, 16, 9 };
                                       // their counter-parts are the
                                       // same +50
      for (unsigned int i=0; i<8; ++i)
        {
                                           // generate a face
          const CellData<2> face = 
            { { vertices[i],      vertices[i+1] ,
                vertices[i+1]+50, vertices[i]+50 }, 2 };
          sub_cell_data.boundary_quads.push_back (face);

                                           // same for the faces
          const CellData<1> edges[4] =
            { { { vertices[i],      vertices[i+1]    }, 2 },
              { { vertices[i]+50,   vertices[i+1]+50 }, 2 },
              { { vertices[i],      vertices[i]+50   }, 2 },
              { { vertices[i+1],    vertices[i+1]+50 }, 2 } };
          for (unsigned int j=0; j<4; ++j)
            sub_cell_data.boundary_lines.push_back (edges[j]);
        }
    }  


                                     // finally top the building
                                     // with four closing cells and
                                     // the vertex at the top
    {
      vertices.push_back (Point<3> (0,0,PhantomGeometry::r1));

      const unsigned int connectivity[4][8]
        = { { 59, 60, 61, 67,   51, 52, 53, 50 },
            { 61, 62, 63, 67,   53, 54, 55, 50 },
            { 67, 63, 64, 65,   50, 55, 56, 57 },
            { 59, 67, 65, 66,   51, 50, 57, 58 }};
      
      for (unsigned int i=0; i<4; ++i)
        {
          CellData<3> cell;
          for (unsigned int j=0; j<8; ++j)
            cell.vertices[j]   = connectivity[i][j];
          cell.material_id   = 0;
          cells.push_back (cell);
        }

                                       // generate boundary
                                       // information for these cells,
                                       // too
      for (unsigned int i=0; i<4; ++i)
        {
          const CellData<2> face = 
            { { connectivity[i][0], connectivity[i][1],
                connectivity[i][2], connectivity[i][3] }, 2 };
          sub_cell_data.boundary_quads.push_back (face);

          const CellData<1> edges[4] =
            { { { connectivity[i][0], connectivity[i][1] }, 2 },
              { { connectivity[i][1], connectivity[i][2] }, 2 },
              { { connectivity[i][2], connectivity[i][3] }, 2 },
              { { connectivity[i][3], connectivity[i][0] }, 2 } };
          for (unsigned int j=0; j<4; ++j)
            sub_cell_data.boundary_lines.push_back (edges[j]);
        }
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

  
  
