// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Create a mesh and attach a function to the signal that is triggered
// whenever we refine the mesh. count how often it is called on each
// processor. (note: this is more than once per refinement cycle because we
// have to re-balance and rebuild the mesh.)

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



int counter = 0;
void
listener()
{
  ++counter;
}


void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<2> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr);
  tr.signals.post_refinement.connect(&listener);

  // try some global refinement
  counter = 0;
  tr.refine_global(1);
  if (myid == 0)
    deallog << "refine_global(1) results in a total of " << counter
            << std::endl;

  counter = 0;
  tr.refine_global(3);
  if (myid == 0)
    deallog << "refine_global(3) results in a total of " << counter
            << std::endl;


  // now also find the bottom left corner of the domain and, on the processor
  // that owns this cell, refine it 4 times
  for (unsigned int i = 0; i < 4; ++i)
    {
      counter                              = 0;
      Triangulation<2>::cell_iterator cell = tr.begin(0);
      while (cell->has_children())
        cell = cell->child(0);
      if (cell->is_locally_owned())
        cell->set_refine_flag();
      tr.execute_coarsening_and_refinement();
      if (myid == 0)
        deallog << "local refinement results in a total of " << counter
                << std::endl;
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      test();
    }
  else
    test();
}
