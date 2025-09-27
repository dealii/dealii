// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test based on parallel_vector_03a but using import function. This test was
// not working because we would only set ghost elements and no local elements
// which resulted in an error during the import

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"

void
check(const unsigned int                                myid,
      const LinearAlgebra::distributed::Vector<double> &v)
{
  v.print(std::cout);
  if (myid == 0)
    {
      AssertThrow(v(10) == 10.0, ExcInternalError());
      AssertThrow(v(11) == 0., ExcInternalError());
      AssertThrow(v(12) == 0., ExcInternalError());
      AssertThrow(v(14) == 14., ExcInternalError());

      AssertThrow(v(5) == 55., ExcInternalError());
    }
  else
    {
      AssertThrow(v(4) == 0., ExcInternalError());
      AssertThrow(v(5) == 55., ExcInternalError());
      AssertThrow(v(6) == 66., ExcInternalError());
    }
}


void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  Assert(numproc == 2, ExcNotImplemented());

  const unsigned int size = 20;
  IndexSet           local_owned(size);
  IndexSet           local_nonzero(size);
  IndexSet           local_relevant(size);
  if (myid == 0)
    {
      local_owned.add_range(0, 10);
      local_nonzero.add_range(5, 10);
      local_relevant = local_owned;
      local_relevant.add_range(10, 13);
      local_relevant.add_range(14, 15);
    }
  else
    {
      local_owned.add_range(10, size);
      local_nonzero.add_range(10, 11);
      local_nonzero.add_range(13, 15);
      local_relevant = local_owned;
      local_relevant.add_range(4, 7);
    }

  LinearAlgebra::distributed::Vector<double> v(local_owned,
                                               local_relevant,
                                               MPI_COMM_WORLD);
  v = 0.;

  // set local values
  for (unsigned int i = 0; i < local_nonzero.n_elements(); ++i)
    v(local_nonzero.nth_index_in_set(i)) = local_nonzero.nth_index_in_set(i);

  // set value from processor which does not own it:
  v(5) = 55.;
  v.compress(VectorOperation::insert);

  // add to value from processor which has it as a ghost
  IndexSet                               ghost_set(size);
  LinearAlgebra::ReadWriteVector<double> rw_vector;
  if (myid == 1)
    {
      ghost_set.add_index(6);
      rw_vector.reinit(ghost_set);
      rw_vector(6) = 60;
    }
  else
    {
      rw_vector.reinit(ghost_set);
    }
  v.import_elements(rw_vector, VectorOperation::add); // 60 + 6
  // compress(insert) used to leave ghosts un-touched which resulted in
  // the wrong 55+55 for this compress(add) operation.

  v.update_ghost_values();

  check(myid, v);

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
