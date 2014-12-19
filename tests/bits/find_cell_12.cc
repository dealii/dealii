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



// find_active_cell_around_point goes into an endless loop, reported
// on the mailing list by Giorgos Kourakos (2014-04-10).

#include "../tests.h"

#include <stdio.h>
#include <cstdlib>

#include <base/quadrature_lib.h>
#include <fe/mapping_q.h>
#include <base/function.h>
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_tools.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <grid/grid_in.h>


#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <sstream>
#include <time.h>

using namespace dealii;

void test()
{
  Triangulation<2> triangulation;

  Point<2> left_bottom(0,-270);
  Point<2> right_top(5000,30);
  std::vector<unsigned int> n_cells;
  n_cells.push_back(10);
  n_cells.push_back(2);


  GridGenerator::subdivided_hyper_rectangle(triangulation,
					    n_cells,
					    left_bottom,
					    right_top,
					    true);

  typename Triangulation<2>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
  for (; cell!=endc; ++cell){
    Point<2> cell_center = cell->center();
    if (abs(cell_center(0) - 1500) < 550){
      cell->set_refine_flag ();
    }
  }

  triangulation.execute_coarsening_and_refinement ();

  Point<2> test_point(250, 195);
  std::cout << "Checking Point " << test_point << std::endl;
  try
    {
      std::pair<typename Triangulation<2>::active_cell_iterator, Point<2> > current_cell =
	GridTools::find_active_cell_around_point(MappingQ1<2>(), triangulation, test_point);

      deallog << "cell: index = " << current_cell.first->index()
	      << " level = " << current_cell.first->level() << std::endl;
      deallog << " pos: " << current_cell.second << std::endl;
    }
  catch (GridTools::ExcPointNotFound<2> &e)
    {
      deallog << "outside" << std::endl;
    }
  deallog << "done" << std::endl;
}

int main (int argc, char **argv)
{
  initlog();

  test();

  return 0;
}
