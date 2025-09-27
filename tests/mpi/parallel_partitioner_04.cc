// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check global_to_local and local_to_global on test case from
// parallel_partitioner_01.cc

#include <deal.II/base/index_set.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/utilities.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;

  const unsigned int set = 200;
  AssertIndexRange(numproc, set - 2);
  const unsigned int local_size  = set - myid;
  unsigned int       global_size = 0;
  unsigned int       my_start    = 0;
  for (unsigned int i = 0; i < numproc; ++i)
    {
      global_size += set - i;
      if (i < myid)
        my_start += set - i;
    }
  // each processor owns some indices and all
  // are ghosting elements from three
  // processors (the second). some entries
  // are right around the border between two
  // processors
  IndexSet local_owned(global_size);
  local_owned.add_range(my_start, my_start + local_size);
  IndexSet local_relevant(global_size);
  local_relevant                 = local_owned;
  unsigned int ghost_indices[10] = {1,
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

  Utilities::MPI::Partitioner v(local_owned, local_relevant, MPI_COMM_WORLD);

  // check locally owned range
  for (unsigned int i = my_start; i < my_start + local_size; ++i)
    {
      AssertDimension(v.global_to_local(i), i - my_start);
      AssertDimension(v.local_to_global(i - my_start), i);
    }

  // check ghost indices
  for (unsigned int i = 0, count = 0; i < 10; ++i)
    if (ghost_indices[i] < my_start ||
        ghost_indices[i] >= my_start + local_size)
      {
        AssertDimension(local_size + count,
                        v.global_to_local(ghost_indices[i]));
        AssertDimension(ghost_indices[i],
                        v.local_to_global(local_size + count));
        ++count;
      }

  // check that loc->glob and glob->loc form an
  // identity operation
  for (unsigned int i = 0; i < local_size + v.n_ghost_indices(); ++i)
    AssertDimension(i, v.global_to_local(v.local_to_global(i)));

  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}
