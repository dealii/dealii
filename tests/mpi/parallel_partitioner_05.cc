// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check the parallel partitioner with an additional index set to describe a
// larger set for which the partitioner provides us with some renumbering
#include <deal.II/base/index_set.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/utilities.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

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
  IndexSet local_relevant(global_size);
  local_relevant                            = local_owned;
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
  local_relevant.add_indices(&ghost_indices[0], &ghost_indices[0] + 10);
  types::global_dof_index before_start = myid > 0 ? my_start - set / 4 : 0;
  types::global_dof_index after_end    = myid < numproc - 1 ?
                                           my_start + local_size + set / 3 :
                                           my_start + local_size;
  if (before_start < my_start)
    local_relevant.add_range(before_start, my_start);
  if (after_end > my_start + local_size)
    local_relevant.add_range(my_start + local_size, after_end);

  Utilities::MPI::Partitioner v(local_owned, local_relevant, MPI_COMM_WORLD);

  std::vector<types::global_dof_index> restricted_indices;
  for (types::global_dof_index i = before_start; i < my_start; ++i)
    if (i % 2 == 0)
      restricted_indices.push_back(i);
  for (types::global_dof_index i = my_start + local_size; i < after_end; ++i)
    if (((i / 4) % 3) == 0)
      restricted_indices.push_back(i);
  restricted_indices.push_back(13);
  restricted_indices.push_back(2 * set);
  restricted_indices.push_back(2 * set + 1);
  IndexSet restricted_set(global_size);
  restricted_set.add_indices(restricted_indices.begin(),
                             restricted_indices.end());
  Utilities::MPI::Partitioner w(local_owned, MPI_COMM_WORLD);
  w.set_ghost_indices(restricted_set, v.ghost_indices());

  IndexSet restricted_set2(global_size);
  restricted_set2.add_index(2);
  if (before_start < my_start)
    restricted_set2.add_range(before_start, my_start);
  Utilities::MPI::Partitioner x(local_owned, MPI_COMM_WORLD);
  x.set_ghost_indices(restricted_set2, v.ghost_indices());

  // print the additional ghost index number to the log stream
  deallog << "Ghost subset in " << v.n_ghost_indices() << " indices: ";
  for (unsigned int i = 0; i < v.ghost_indices_within_larger_ghost_set().size();
       ++i)
    deallog << '[' << v.ghost_indices_within_larger_ghost_set()[i].first << ", "
            << v.ghost_indices_within_larger_ghost_set()[i].second << ") ";
  deallog << std::endl;

  deallog << "Ghost subset in " << w.n_ghost_indices() << " indices: ";
  for (unsigned int i = 0; i < w.ghost_indices_within_larger_ghost_set().size();
       ++i)
    deallog << '[' << w.ghost_indices_within_larger_ghost_set()[i].first << ", "
            << w.ghost_indices_within_larger_ghost_set()[i].second << ") ";
  deallog << std::endl;

  deallog << "Ghost subset in " << x.n_ghost_indices() << " indices: ";
  for (unsigned int i = 0; i < x.ghost_indices_within_larger_ghost_set().size();
       ++i)
    deallog << '[' << x.ghost_indices_within_larger_ghost_set()[i].first << ", "
            << x.ghost_indices_within_larger_ghost_set()[i].second << ") ";
  deallog << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv);
  MPILogInitAll                    log;
  test();
}
