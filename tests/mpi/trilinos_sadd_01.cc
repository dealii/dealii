// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check correct behavior of sadd of Trilinos vectors
// if they have different Epetra maps

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_vector.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "../tests.h"

void
test()
{
  const int n_proc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const int my_id  = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  // All processes should own 10 entries
  const int entries_per_process = 10;

  IndexSet  locally_owned(entries_per_process * n_proc);
  const int begin_index = my_id * entries_per_process;
  const int end_index   = (my_id + 1) * entries_per_process;
  locally_owned.add_range(begin_index, end_index);

  IndexSet  locally_relevant(entries_per_process * n_proc);
  const int local_begin = std::max(0, begin_index - entries_per_process / 2);
  const int local_end   = entries_per_process * n_proc;
  locally_relevant.add_range(local_begin, local_end);

  TrilinosWrappers::MPI::Vector ghosted, distributed;
  distributed.reinit(locally_owned, MPI_COMM_WORLD);
  ghosted.reinit(locally_owned, locally_relevant, MPI_COMM_WORLD);

  // set the 'distributed' vector to all ones and store its results in
  // 'ghosted'
  distributed = 1.;
  ghosted     = distributed;

  // then multiply 'distributed' by two and add 'ghosted' to it again
  distributed *= 2;
  distributed.add(ghosted, true);

  // assign the result, which should contain all 3s to 'ghosted'
  ghosted = distributed;

  if (my_id == 0)
    {
      deallog << "Distributed:" << std::endl;
      for (int i = begin_index; i < end_index; ++i)
        deallog << i << ": " << distributed(i) << std::endl;

      deallog << "Ghosted:" << std::endl;
      for (int i = local_begin; i < local_end; ++i)
        deallog << i << ": " << ghosted(i) << std::endl;
    }

  // verify correct value
  for (int i = begin_index; i < end_index; ++i)
    Assert(distributed(i) == 3, ExcInternalError());

  for (int i = local_begin; i < local_end; ++i)
    Assert(ghosted(i) == 3, ExcInternalError());
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
