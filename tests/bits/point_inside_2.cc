//----------------------------  point_inside_2.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003 by the deal.II authors and Anna Schneebeli
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  point_inside_2.cc  ---------------------------


// check TriaAccessor<3>::point_inside for a cell that is _not_
// aligned with the coordinate axes
//
// this program is a modified version of one by Joerg Weimar,
// TU Braunschweig

#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <fstream>


void check ()
{
				   // use a rather complicated way to
				   // generate a new mesh with one
				   // moved vertex -- one could move
				   // it by hand in the original
				   // triangulation, but this is how
				   // the original testcase came in,
				   // and it tests some more stuff, so
				   // why change things
  Triangulation<3> triangulation;
  
				   // we generate a hyper_cube and
				   // modify one vertex.
  GridGenerator::hyper_cube (triangulation);
  
				   // Now get the cell
  Triangulation<3>::cell_iterator cell = triangulation.begin();
    
  std::vector<CellData<3> > cells (1, CellData<3>());
				   // get the indices
  for (unsigned int j=0; j<8; ++j) {       
    cells[0].vertices[j] = cell->vertex_index(j);
  }
    
				   // get the vertex coordinates
  Point<3> vertices[8];
  for (unsigned int j=0; j<8; ++j) {       
    vertices[j] = cell->vertex(j);
  }
				   // modify one vertex.
  vertices[0] = Point<3>(-1,0,0);
    
  SubCellData boundary_info;
  Triangulation<3> triangulation2;
  
				   // create a new triangulation.
  triangulation2.create_triangulation (
    std::vector<Point<3> >(&vertices[0],&vertices[8]),
    cells , boundary_info
  );

  Triangulation<3>::cell_iterator cell2 = triangulation2.begin();
                                      
				   // and test it.                                   
  double testcoord[14][3] = {{0.5,0.5,0.5},
			     {2,0.5,0.5},
			     {0.5,2,0.5},
			     {0.5,0.5,2},
			     {-2,0.5,0.5},
			     {0.5,-2,0.5},
			     {0.5,0.5,-2},
			     {0.9,0.9,0.9},
                            {1.0,0.5,0.5},
                            {0.9999999,0.5,0.5},
                            {1.0000001,0.5,0.5},
                            {-0.1,0.1,0.1},
                            {-0.24,0.5,0.5},
                            {-0.26,0.5,0.5}
                            };
    
    int expected[] = {1,0,0,0,0,0,0,1,1,1,0,1,1,0};
    for (int i=0; i<14; i++)
      {
        const Point<3> testpoint(testcoord[i][0],
				 testcoord[i][1],
				 testcoord[i][2]);
        bool res = cell2->point_inside(testpoint);
        deallog << testpoint << " \t inside " << res 
                << " expected " << expected[i] << std::endl;
        Assert (res == expected[i], ExcInternalError());
      }    
}


int main () 
{
  std::ofstream logfile("point_inside_2.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  check ();
}
