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


// find_active_cell_around_point() stops working if some cells are
// refined in 3d.  this is caused by a bug in
// GridTools::find_cells_adjacent_to_vertex that got fixed in r 25562.
//
// the bug here is the same as in find_cell_6 but when calling the
// function with hp:: arguments

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>

#include <deal.II/fe/mapping_q1.h>

bool inside(Triangulation<3> &tria, Point<3> &p)
{

  for (Triangulation<3>::cell_iterator cell = tria.begin(0);
       cell != tria.end(0); ++cell)
    if ( cell->point_inside (p) )
      return true;

  return false;
}

void check2 ()
{
  Triangulation<3> tria;
  GridIn<3> gridIn;
  gridIn.attach_triangulation (tria);
  std::ifstream inputFile(SOURCE_DIR "/grids/grid.inp");
  gridIn.read_ucd (inputFile);

  Point<3> p2(304.767,-57.0113,254.766);

  int idx=0;
  for (Triangulation<3>::active_cell_iterator  cell = tria.begin_active();  cell!=tria.end();  ++cell, ++idx)
    {
      if (idx==21)
        cell->set_refine_flag();
    }
  tria.execute_coarsening_and_refinement ();

  deallog << inside(tria, p2) << std::endl;

  hp::MappingCollection<3> mappings;
  mappings.push_back (MappingQ1<3>());
  mappings.push_back (MappingQ1<3>());

  hp::FECollection<3> fes;
  fes.push_back (FE_Q<3>(1));
  fes.push_back (FE_Q<3>(1));

  hp::DoFHandler<3> dof_handler (tria);
  dof_handler.distribute_dofs (fes);

  GridTools::find_active_cell_around_point(mappings, dof_handler, p2); //triggered exception
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      check2();
    }
  catch (const std::exception &exc)
    {
      // we shouldn't get here...
      deallog << "Caught an error..." << std::endl;
      deallog << exc.what() << std::endl;
    }
}



