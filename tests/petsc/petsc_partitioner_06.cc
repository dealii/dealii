// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for the Partitioner with a smaller ghost index set within a larger one
// regarding the export_to_ghosted_array() calls
// Copy-pasted from tests/mpi/parallel_partitioner_06.c
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/petsc_communication_pattern.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  Assert(numproc > 2, ExcNotImplemented());

  const unsigned int set = 200;
  AssertIndexRange(numproc, set - 2);
  const unsigned int      local_size  = set - myid;
  types::global_dof_index global_size = 0;
  types::global_dof_index my_start    = 0;
  for (unsigned int i = 0; i < numproc; ++i)
    {
      global_size += set - i;
      if (i < myid)
        my_start += set - i;
    }

  // each processor owns some indices and all are ghosting elements from three
  // processors (the second). some entries are right around the border between
  // two processors
  IndexSet local_owned(global_size);
  local_owned.add_range(my_start, my_start + local_size);
  IndexSet local_relevant_1(global_size), local_relevant_2(global_size);
  local_relevant_1                          = local_owned;
  types::global_dof_index ghost_indices[10] = {1,
                                               2,
                                               13,
                                               set - 2,
                                               set - 1,
                                               set,
                                               set + 1,
                                               2 * set,
                                               2 * set + 1,
                                               2 * set + 3};
  local_relevant_1.add_indices(&ghost_indices[0], ghost_indices + 10);
  if (myid > 0)
    local_relevant_1.add_range(my_start - 10, my_start);
  if (myid < numproc - 1)
    local_relevant_1.add_range(my_start + local_size,
                               my_start + local_size + 10);

  local_relevant_2 = local_owned;
  local_relevant_2.add_indices(&ghost_indices[0], ghost_indices + 10);
  if (myid > 0)
    local_relevant_2.add_index(my_start - 10);
  if (myid < numproc - 1)
    local_relevant_2.add_index(my_start + local_size + 9);

  IndexSet local_relevant_3(global_size);
  local_relevant_3.add_index(2);
  if (myid > 0 && my_start > 0)
    local_relevant_3.add_range(my_start - 10, my_start);

  PETScWrappers::Partitioner v, w, x;
  v.reinit(local_owned, local_relevant_1, MPI_COMM_WORLD);
  w.reinit(local_owned, local_relevant_2, v.ghost_indices(), MPI_COMM_WORLD);
  x.reinit(local_owned, local_relevant_3, v.ghost_indices(), MPI_COMM_WORLD);

  // set up a locally owned array with some entries
  std::vector<unsigned int> locally_owned_data_field(local_size);
  for (unsigned int i = 0; i < local_size; ++i)
    locally_owned_data_field[i] = my_start + i;
  const std::vector<unsigned int> &locally_owned_data(locally_owned_data_field);

  // set up a ghost array
  std::vector<unsigned int> ghosts(v.n_ghost_indices());

  // send the full array
  v.export_to_ghosted_array_start(make_array_view(locally_owned_data),
                                  make_array_view(ghosts));
  v.export_to_ghosted_array_finish(make_array_view(locally_owned_data),
                                   make_array_view(ghosts));
  deallog << "All ghosts: ";
  for (unsigned int i = 0; i < ghosts.size(); ++i)
    deallog << ghosts[i] << ' ';
  deallog << std::endl;

  // send only the array in w
  std::fill(ghosts.begin(), ghosts.end(), 0);

  w.export_to_ghosted_array_start(make_array_view(locally_owned_data),
                                  make_array_view(ghosts));


  // start a second send operation for the x partitioner in parallel to make
  // sure communication does not get messed up
  std::vector<unsigned int> ghosts2(x.n_ghost_indices());

  x.export_to_ghosted_array_start(make_array_view(locally_owned_data),
                                  make_array_view(ghosts2));

  w.export_to_ghosted_array_finish(make_array_view(locally_owned_data),
                                   make_array_view(ghosts));
  deallog << "Ghosts on reduced 1: ";
  for (unsigned int i = 0; i < ghosts.size(); ++i)
    deallog << ghosts[i] << ' ';
  deallog << std::endl;

  std::fill(ghosts.begin(), ghosts.end(), 0);

  x.export_to_ghosted_array_start(make_array_view(locally_owned_data),
                                  make_array_view(ghosts));
  x.export_to_ghosted_array_finish(make_array_view(locally_owned_data),
                                   make_array_view(ghosts));
  deallog << "Ghosts on reduced 2: ";
  for (unsigned int i = 0; i < ghosts.size(); ++i)
    deallog << ghosts[i] << ' ';
  deallog << std::endl;

  x.export_to_ghosted_array_finish(make_array_view(locally_owned_data),
                                   make_array_view(ghosts2));
  deallog << "Ghosts on reduced 2 without excess entries: ";
  for (unsigned int i = 0; i < ghosts2.size(); ++i)
    deallog << ghosts2[i] << ' ';
  deallog << std::endl;

  x.export_to_ghosted_array_start(make_array_view(locally_owned_data),
                                  make_array_view(ghosts));
  x.export_to_ghosted_array_finish(make_array_view(locally_owned_data),
                                   make_array_view(ghosts));
  deallog << "Ghosts on reduced 2: ";
  for (unsigned int i = 0; i < ghosts.size(); ++i)
    deallog << ghosts[i] << ' ';
  deallog << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv);
  MPILogInitAll                    log;
  test();
}
