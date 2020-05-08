// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// similar to parallel_sparse_vector_03.cc, but make sure
// compress(insert) zeroes out ghosts in Release mode

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "../tests.h"

void
check(const unsigned int                                myid,
      const LinearAlgebra::distributed::Vector<double> &v)
{
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
  for (unsigned int i = 0; i < local_nonzero.n_elements(); i++)
    v(local_nonzero.nth_index_in_set(i)) = local_nonzero.nth_index_in_set(i);

  // set value from processor which does not own it:
  v(5) = 55.;
  v.compress(VectorOperation::insert);

  // add to value from processor which has it as a ghost
  if (myid == 1)
    v(6) = 60;
  v.compress(VectorOperation::add); // 60 + 6
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
