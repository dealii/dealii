// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// This test provides a hp::DoFHandler with two FE_Q(2) elements
// assigned on a p::d::Triangulation consisting of 64x64 cells, on which
// active FE indices are mostly randomly distributed. This DoF
// distribution test repeats for an increasing amount of processors.
//
// At some point, this test failed to provide a consistent number of
// degrees of freedom for different numbers of processors assigned.


#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include <iostream>

#include "../tests.h"


template <int dim>
void
test(MPI_Comm mpi_communicator)
{
  parallel::distributed::Triangulation<dim> triangulation(mpi_communicator);

  hp::FECollection<dim> fe;
  fe.push_back(FE_Q<dim>(2));
  fe.push_back(FE_Q<dim>(2));

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(4);

  DoFHandler<dim> dof_handler(triangulation);

  // set active_fe_index mostly randomly
  for (const auto &cell : dof_handler.active_cell_iterators() |
                            IteratorFilters::LocallyOwnedCell())
    cell->set_active_fe_index(cell->active_cell_index() % fe.size());

  dof_handler.distribute_dofs(fe);

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    deallog << "n_procs/n_dofs: "
            << Utilities::MPI::n_mpi_processes(mpi_communicator) << '/'
            << dof_handler.n_dofs() << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));
  if (myid == 0)
    initlog();

  // Create communicators with larger and larger subsets of
  // processors, and then see what happens with the number of DoFs we
  // get on the (always same) mesh created in the test() function
  for (unsigned int n = 1; n <= Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
       ++n)
    {
      // Colorize the current processor: zero if we want to include it
      // in the communication, one otherwise
      const int color =
        (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) < n ? 0 : 1);

      // Then split the communicator into two, based on the two
      // colors. We'll only do something with that communicator that
      MPI_Comm subset_comm;
      MPI_Comm_split(MPI_COMM_WORLD, color, /*key=*/0, &subset_comm);

      if (color == 0)
        test<2>(subset_comm);

      // destroy the communicators previously created
      MPI_Comm_free(&subset_comm);
    }
}
