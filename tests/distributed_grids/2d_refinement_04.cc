// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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



// Test interaction with p4est with a simple grid in 2d. here, we
// check that if we refine a square once, and then one of the children
// once again, that we get 7 cells
//
// at the time of writing this test, the results for this testcase
// were erratic and apparently non-deterministic. the actual cause was
// an uninitialized variable, fixed in revision 16414.

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

#include <fstream>
#include <cstdlib>


template<int dim>
void test(std::ostream & /*out*/)
{
  for (unsigned int i=0; i<GeometryInfo<dim>::max_children_per_cell; ++i)
    {
      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::hyper_cube(tr);

      deallog << i << ' ' << tr.n_active_cells()
              << std::endl;

      tr.refine_global (1);

      deallog << i << ' ' << tr.n_active_cells()
              << std::endl;

      Assert (tr.n_active_cells() == 4,
              ExcInternalError());

      typename Triangulation<dim>::active_cell_iterator
      cell = tr.begin_active();
      std::advance (cell, i);
      cell->set_refine_flag();
      tr.execute_coarsening_and_refinement ();

      deallog << i << ' ' << tr.n_active_cells()
              << std::endl;

      Assert (tr.n_active_cells() == 7,
              ExcInternalError());
    }
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("2d");
  test<2>(logfile);
  deallog.pop();


}
