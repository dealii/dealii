// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test AlignedVector::replicate_across_communicator().
//
// Check what happens if we call this function on an empty object.

#include <deal.II/base/aligned_vector.h>

#include "../tests.h"


void
test()
{
  const MPI_Comm     communicator = MPI_COMM_WORLD;
  const unsigned int root         = 1;
  Assert(root < Utilities::MPI::n_mpi_processes(communicator),
         ExcInternalError());

  // Create an empty object and replicate it.
  AlignedVector<int> avec;
  avec.replicate_across_communicator(communicator, root);

  deallog << "On process " << Utilities::MPI::this_mpi_process(communicator)
          << ": " << avec.size() << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test();
}
