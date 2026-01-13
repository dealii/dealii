// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test reading from ghost entries of vectors

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

  fully_distributed_vector = 1.0;

  // set local values
  fully_distributed_vector(myid * 2)     = myid * 2.0;
  fully_distributed_vector(myid * 2 + 1) = myid * 2.0 + 1.0;
  fully_distributed_vector.compress(VectorOperation::insert);
  fully_distributed_vector *= 2.0;

  typename LA::MPI::Vector vector_with_ghost_entries(local_active,
                                                     local_relevant,
                                                     MPI_COMM_WORLD);
  vector_with_ghost_entries = fully_distributed_vector;
  Assert(vector_with_ghost_entries.has_ghost_elements(), ExcInternalError());

  // Each process has entry 1 as a ghost (or actually owns it). Check this:
  deallog << "ghosted value: "
          << get_real_assert_zero_imag(vector_with_ghost_entries(1))
          << " (should be 2.0)" << std::endl;

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
