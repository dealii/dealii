// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check AffineConstraints<double>.distribute() for a petsc vector
//
// like _01, but for a block vector. this has the additional complication that
// (at a global level) the set of indices owned by this processor is not
// contiguous

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/petsc_block_vector.h>

#include "../tests.h"



void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_processes =
    Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  // create a vector that consists of elements indexed from 0 to n
  PETScWrappers::MPI::BlockVector vec(2,
                                      MPI_COMM_WORLD,
                                      100 * n_processes,
                                      100);
  vec.block(0).reinit(MPI_COMM_WORLD, 100 * n_processes, 100);
  vec.block(1).reinit(MPI_COMM_WORLD, 100 * n_processes, 100);
  vec.collect_sizes();
  AssertThrow(vec.block(0).locally_owned_size() == 100, ExcInternalError());
  AssertThrow(vec.block(0).local_range().first == 100 * myid,
              ExcInternalError());
  AssertThrow(vec.block(0).local_range().second == 100 * myid + 100,
              ExcInternalError());
  AssertThrow(vec.block(1).locally_owned_size() == 100, ExcInternalError());
  AssertThrow(vec.block(1).local_range().first == 100 * myid,
              ExcInternalError());
  AssertThrow(vec.block(1).local_range().second == 100 * myid + 100,
              ExcInternalError());

  IndexSet locally_owned(vec.size());
  locally_owned.add_range(100 * myid, 100 * myid + 100);
  locally_owned.add_range(vec.block(0).size() + 100 * myid,
                          vec.block(0).size() + 100 * myid + 100);
  AssertThrow(vec.locally_owned_elements() == locally_owned,
              ExcInternalError());

  if (myid == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      test();
    }
  else
    test();
}
