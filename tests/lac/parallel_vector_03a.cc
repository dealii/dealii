// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// similar to parallel_sparse_vector_03.cc, but make sure
// compress(insert) zeroes out ghosts in Release mode

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "../tests.h"

void
check(const unsigned int myid, const LinearAlgebra::ReadWriteVector<double> &v)
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

  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> v(
    local_owned, local_relevant, MPI_COMM_WORLD);
  v = 0.;

  // set local values
  IndexSet indexset_1(local_owned);
  indexset_1.add_index(5);
  indexset_1.compress();
  LinearAlgebra::ReadWriteVector<double> rw_vector(indexset_1);
  for (unsigned int i = 0; i < local_nonzero.n_elements(); ++i)
    rw_vector(local_nonzero.nth_index_in_set(i)) =
      local_nonzero.nth_index_in_set(i);

  // set value from processor which does not own it:
  rw_vector(5) = 55.;
  v.import_elements(rw_vector, VectorOperation::insert);

  // add to value from processor which has it as a ghost
  // Because of limitation in import, the ReadWriteVector needs to have locally
  // owned elements
  IndexSet                               indexset_2(local_owned.size());
  LinearAlgebra::ReadWriteVector<double> rw_add;
  if (myid == 1)
    {
      indexset_2.add_index(6);
      rw_add.reinit(indexset_2);
      rw_add(6) = 60;
    }
  else
    {
      rw_add.reinit(indexset_2);
    }
  v.import_elements(rw_add, VectorOperation::add); // 60 + 6
  // compress(insert) used to leave ghosts un-touched which resulted in
  // the wrong 55+55 for this compress(add) operation.

  v.update_ghost_values();

  IndexSet indexset_3(size);
  if (myid == 0)
    {
      indexset_3.add_index(10);
      indexset_3.add_index(11);
      indexset_3.add_index(12);
      indexset_3.add_index(14);

      indexset_3.add_index(5);
    }
  else
    {
      indexset_3.add_index(4);
      indexset_3.add_index(5);
      indexset_3.add_index(6);
    }

  rw_vector.reinit(indexset_3);
  rw_vector.import_elements(v, VectorOperation::insert);

  check(myid, rw_vector);

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
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}
