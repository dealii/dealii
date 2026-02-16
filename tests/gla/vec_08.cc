// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// We allow zeroing out vectors that have ghost elements, as the only
// operation that modifies vectors that have ghost entries. Check that
// that is possible, and that the results are correct for both locally
// owned and ghost entries.

#include <deal.II/base/index_set.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/generic_linear_algebra.h>

#include <iostream>
#include <vector>

#include "../tests.h"

#include "gla.h"

template <class LA>
void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;

  // each processor owns 2 indices and all
  // are ghosting Element 1 (the second)

  IndexSet local_active(numproc * 2);
  local_active.add_range(myid * 2,
                         myid * 2 + 2); // Process k owns entries 2k, 2k+1

  IndexSet local_relevant(numproc * 2);
  local_relevant.add_range(1, 2); // All processes also have entry 1 as ghost

  typename LA::MPI::Vector fully_distributed_vector(local_active,
                                                    MPI_COMM_WORLD);

  // Create a vector with nonzero entries for all locally owned entries:
  fully_distributed_vector(myid * 2)     = myid * 2.0 + 1.0;
  fully_distributed_vector(myid * 2 + 1) = myid * 2.0 + 2.0;
  fully_distributed_vector.compress(VectorOperation::insert);

  // Then create a vector with ghost entries that matches the fully
  // distributed one:
  typename LA::MPI::Vector vector_with_ghost_entries(local_active,
                                                     local_relevant,
                                                     MPI_COMM_WORLD);
  vector_with_ghost_entries = fully_distributed_vector;
  Assert(vector_with_ghost_entries.has_ghost_elements(), ExcInternalError());

  Assert(vector_with_ghost_entries(myid * 2) != 0, ExcInternalError());
  Assert(vector_with_ghost_entries(myid * 2 + 1) != 0, ExcInternalError());
  Assert(vector_with_ghost_entries(1) != 0, ExcInternalError());

  // Now zero out the whole vector. The result needs to be zero in
  // both the locally owned entries (2k, 2k+1) and in the ghost entry
  // (entry 1):
  vector_with_ghost_entries = 0;
  Assert(vector_with_ghost_entries(myid * 2) == 0, ExcInternalError());
  Assert(vector_with_ghost_entries(myid * 2 + 1) == 0, ExcInternalError());
  Assert(vector_with_ghost_entries(1) == 0, ExcInternalError());

  // done
  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  {
    deallog.push("PETSc");
    test<LA_PETSc>();
    deallog.pop();
    deallog.push("Trilinos");
    test<LA_Trilinos>();
    deallog.pop();
  }
}
