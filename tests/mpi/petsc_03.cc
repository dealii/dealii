// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// PETScWrappers: document bug when using GrowingVectorMemory
// at the end of a run.

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/vector_memory.h>

#include "../tests.h"


template <class v>
void
testit()
{
  GrowingVectorMemory<v> m;
  m.free(m.alloc());
}



void
test()
{
  unsigned int myid     = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "Running on " << numprocs << " CPU(s)." << std::endl;

  testit<PETScWrappers::MPI::Vector>();
  testit<PETScWrappers::MPI::BlockVector>();

  if (myid == 0)
    deallog << "done" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();

      test();
    }
  else
    test();
}
