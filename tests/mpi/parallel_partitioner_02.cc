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


// check n_import_indices on test case from parallel_partitioner_01

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
  local_relevant = local_owned;
  if (myid < 6)
    {
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
    }

  Utilities::MPI::Partitioner v(local_owned, local_relevant, MPI_COMM_WORLD);

  // check number of import indices everywhere
  // (counted the above) times the number of
  // processors which have these as ghosts
  const unsigned int n_procs_with_ghosts = std::min(numproc - 1, 5U);
  if (myid == 0)
    {
      AssertDimension(v.n_import_indices(), 5 * n_procs_with_ghosts);
    }
  else if (myid == 1)
    {
      AssertDimension(v.n_import_indices(), 2 * n_procs_with_ghosts);
    }
  else if (myid == 2)
    {
      AssertDimension(v.n_import_indices(), 3 * n_procs_with_ghosts);
    }
  else
    {
      AssertDimension(v.n_import_indices(), 0);
    }


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
