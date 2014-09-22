// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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


// eliminate_unrefined_islands has this nasty effect that one marked
// cell can trigger a second cell that had been visited previously in
// the loop over all cells to be marked, etc. in
// prepare_coarsen_and_refinement, this meant that we may have to
// iterate the outer loop quite a number of times, which is not
// efficient.
//
// this should now have become better, since upon marking one cell we
// also investigate all its neighbors as well


char logname[] = "output";


#include "../tests.h"


#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <fstream>
#include <iostream>


using namespace dealii;



void test ()
{
  // generate a 100x3 mesh
  Triangulation<2> triangulation (Triangulation<2>::eliminate_unrefined_islands);
  std::vector<unsigned int> ref(2);
  ref[0] = 100;
  ref[1] = 3;
  GridGenerator::subdivided_hyper_rectangle (triangulation, ref,
                                             Point<2>(), Point<2>(100,3));

  // refine all cells at the lower
  // boundary. we then have 600 cells
  for (Triangulation<2>::cell_iterator
       cell = triangulation.begin();
       cell != triangulation.end(); ++cell)
    if (cell->center()[1] < 1)
      cell->set_refine_flag ();
  triangulation.execute_coarsening_and_refinement ();
  deallog << "n_active_cells = " << triangulation.n_active_cells()
          << std::endl;

  // now mark all cells at the top
  // boundary for refinement with the
  // exception of the top right
  // one. this means that the only
  // cell that qualifies the island
  // condition is the one in the
  // center at the left boundary
  // (bottom neighbor already
  // refined, top neighbor to be
  // refined). but upon it being
  // flagged, its right neighbor will
  // also qualify, and then the right
  // neighbor of that one, etc. in
  // total, the old algorithm needs
  // 102 iterations through the
  // prepare_c_and_r loop (one for
  // each of the cells in the middle
  // stripe, plus two for apparently
  // other things). after the changes
  // to tria.cc, we now need 2
  // iterations
  for (Triangulation<2>::cell_iterator
       cell = triangulation.begin();
       cell != triangulation.end(); ++cell)
    if (cell->center()[1] > 2)
      if (cell->center()[0] < 99)
        cell->set_refine_flag ();
  triangulation.execute_coarsening_and_refinement ();

  // output the new number of
  // cells. should now be 1200 since
  // we have refined every cell in
  // the mesh. unfortunately, there
  // is no way to test actually test
  // the number of iterations in
  // prepare_c_and_r :-(
  deallog << "n_active_cells = " << triangulation.n_active_cells()
          << std::endl;
}



int main ()
{
  std::ofstream logfile(logname);
  logfile.precision (3);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  try
    {
      test();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
