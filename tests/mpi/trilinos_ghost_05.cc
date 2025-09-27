// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check correct behavior of Trilinos ghosted vectors
// check if assignment from a normal vector works correctly and updates the
// ghost values

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_vector.h>

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

  unsigned int ghostel = (numproc > 1) ? 2 : 1;

  // each processor owns 2 indices and all
  // are ghosting one element
  IndexSet local_active(numproc * 2);
  local_active.add_range(myid * 2, myid * 2 + 2);
  IndexSet local_relevant(numproc * 2);
  local_relevant = local_active;
  local_relevant.add_range(ghostel, ghostel + 1);

  TrilinosWrappers::MPI::Vector x(local_active, MPI_COMM_WORLD);
  TrilinosWrappers::MPI::Vector v(local_relevant, MPI_COMM_WORLD);

  // set local values
  x(myid * 2)     = myid * 2.0;
  x(myid * 2 + 1) = myid * 2.0 + 1.0;

  // transfer to ghosted vector v and check
  x.compress(VectorOperation::insert);
  v = x;

  Assert(v(myid * 2) == myid * 2.0, ExcInternalError());
  Assert(v(myid * 2 + 1) == myid * 2.0 + 1.0, ExcInternalError());
  Assert(v(ghostel) == ghostel, ExcInternalError());

  // change x, transfer, and check again
  x *= 2.0;
  v = x;

  Assert(v(myid * 2) == myid * 4.0, ExcInternalError());
  Assert(v(myid * 2 + 1) == myid * 4.0 + 2.0, ExcInternalError());
  Assert(v(ghostel) == 2.0 * ghostel, ExcInternalError());

  if (myid == 0)
    {
      deallog << myid * 2 << ':' << v(myid * 2) << std::endl;
      deallog << myid * 2 + 1 << ':' << v(myid * 2 + 1) << std::endl;
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
