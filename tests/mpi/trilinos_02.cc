// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check Trilinos has_ghost_elements()

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"

void
check(TrilinosWrappers::MPI::Vector &v, bool ghost)
{
  Assert(v.has_ghost_elements() == ghost, ExcMessage("wrong ghost elements"));
}

void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;

  // each processor owns 2 indices and all
  // are ghosting element 1 (the second)
  IndexSet local_active(numproc * 2);
  local_active.add_range(myid * 2, myid * 2 + 2);
  IndexSet local_relevant(numproc * 2);
  local_relevant = local_active;
  local_relevant.add_range(1, 2);

  TrilinosWrappers::MPI::Vector v(local_active, MPI_COMM_WORLD);
  check(v, false);

  TrilinosWrappers::MPI::Vector v2(local_relevant, MPI_COMM_WORLD);
  check(v2, true);

  TrilinosWrappers::MPI::Vector v3 = v;
  check(v3, false);

  v3.reinit(v2);
  check(v3, true);

  v3.reinit(local_active);
  check(v3, false);

  v3.reinit(local_relevant);
  check(v3, true);

  TrilinosWrappers::MPI::Vector v4 = v2;
  check(v4, true);

  v4 = v;
  check(v4, true); // this only copies contents!

  v4 = v2;
  check(v4, true);

  TrilinosWrappers::MPI::Vector v5(v2);
  check(v5, true);

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
