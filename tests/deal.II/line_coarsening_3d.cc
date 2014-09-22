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
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iomanip>
#include <iomanip>


std::ofstream logfile("output");


// create a triangulation of a cylinder, where the line along the axis is common
// to all cells, thus n_cells are neighboring at this line. in delete_children,
// this is not accounted for (for n_cells>5) and the child lines are deleted,
// although they are still needed.
void create_star_structured_cylinder (Triangulation<3> &coarse_grid,
                                      const unsigned int n_cells)
{
  Assert(n_cells>1, ExcNotImplemented());

  std::vector<Point<3> > points(2*(1+2*n_cells));
  points[0] = Point<3>();
  points[1] = Point<3>(1,0,0);
  for (unsigned int i=0; i<2*n_cells-1; ++i)
    {
      points[2+i]=Point<3>(std::cos(numbers::PI/n_cells*(i+1)),
                           std::sin(numbers::PI/n_cells*(i+1)),
                           0);
    }

  for (unsigned int i=0; i<2*n_cells+1; ++i)
    {
      points[1+2*n_cells+i]=points[i]+Point<3>(0,0,-1);
    }

  std::vector<CellData<3> > cells(n_cells);

  for (unsigned int c=0; c<n_cells; ++c)
    {
      cells[c].vertices[0]   = 0;
      cells[c].vertices[1]   = 1+2*c;
      cells[c].vertices[2]   = 2+2*c;
      cells[c].vertices[3]   = (3+2*c)%(2*n_cells);
      cells[c].vertices[4]   = 0+2*n_cells+1;
      cells[c].vertices[5]   = 1+2*c+2*n_cells+1;
      cells[c].vertices[6]   = 2+2*c+2*n_cells+1;
      cells[c].vertices[7]   = (3+2*c)%(2*n_cells)+2*n_cells+1;
    }
  // finally generate a triangulation
  // out of this
  coarse_grid.create_triangulation_compatibility (points, cells, SubCellData());
}


void check ()
{
  const unsigned int dim=3;

  // create tria
  Triangulation<dim> tria;
  create_star_structured_cylinder(tria,6);
  // out of the six cells, refine the first and
  // fourth
  Triangulation<3>::active_cell_iterator cell=tria.begin_active();
  cell->set_refine_flag();
  ++cell;
  ++cell;
  ++cell;
  cell->set_refine_flag();
  tria.execute_coarsening_and_refinement ();
  // write grid to file
  GridOut go;
  go.set_flags (GridOutFlags::Ucd(true));
  go.write_ucd(tria,logfile);
  // coarsen the first cell again
  cell=tria.begin_active(1);
  for (unsigned int c=0; c<8; ++c)
    {
      cell->set_coarsen_flag();
      ++cell;
    }
  tria.execute_coarsening_and_refinement ();

  go.write_ucd(tria,logfile);
}


int main()
{
  deallog.attach(logfile);
  deallog.depth_console(0);
  check();
  return 0;
}

