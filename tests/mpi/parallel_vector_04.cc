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


// check that operator= resets ghosts, both if they have been set and if they
// have not been set

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>

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


  // each processor owns 2 indices and all
  // are ghosting element 1 (the second)
  IndexSet local_owned(numproc * 2);
  local_owned.add_range(myid * 2, myid * 2 + 2);
  IndexSet local_relevant(numproc * 2);
  local_relevant = local_owned;
  local_relevant.add_range(1, 2);

  LinearAlgebra::distributed::Vector<double> v(local_owned,
                                               local_relevant,
                                               MPI_COMM_WORLD);

  // set local values and check them
  v(myid * 2)     = myid * 2.0;
  v(myid * 2 + 1) = myid * 2.0 + 1.0;

  v.compress(VectorOperation::insert);
  v *= 2.0;

  AssertThrow(v(myid * 2) == myid * 4.0, ExcInternalError());
  AssertThrow(v(myid * 2 + 1) == myid * 4.0 + 2.0, ExcInternalError());

  // set ghost dof on remote process, no
  // compress called
  if (myid > 0)
    v(1) = 7;

  AssertThrow(v(myid * 2) == myid * 4.0, ExcInternalError());
  AssertThrow(v(myid * 2 + 1) == myid * 4.0 + 2.0, ExcInternalError());

  if (myid > 0)
    AssertThrow(v(1) == 7.0, ExcInternalError());

  // reset to zero
  v = 0;

  AssertThrow(v(myid * 2) == 0., ExcInternalError());
  AssertThrow(v(myid * 2 + 1) == 0., ExcInternalError());

  // check that everything remains zero also
  // after compress
  v.compress(VectorOperation::add);

  AssertThrow(v(myid * 2) == 0., ExcInternalError());
  AssertThrow(v(myid * 2 + 1) == 0., ExcInternalError());

  // set element 1 on owning process to
  // something nonzero
  if (myid == 0)
    v(1) = 2.;
  if (myid > 0)
    AssertThrow(v(1) == 0., ExcInternalError());

  // check that all processors get the correct
  // value again, and that it is erased by
  // operator=
  v.update_ghost_values();

  AssertThrow(v(1) == 2.0, ExcInternalError());

  v = 0;
  AssertThrow(v(1) == 0.0, ExcInternalError());

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
