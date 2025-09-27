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



// PETScWrappers: document bug with PETSc SparseMatrix. If only one CPU
// does matrix-assembly, it calls compress() inside and the others don't.
// We should implement this like in PETSc::MPI::Vector.

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_sparse_matrix.h>

#include "../tests.h"


void
test()
{
  unsigned int myid     = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "Running on " << numprocs << " CPU(s)." << std::endl;

  DynamicSparsityPattern csp(2);
  csp.add(0, 0);
  csp.add(1, 1);

  PETScWrappers::MPI::SparseMatrix     mat;
  std::vector<types::global_dof_index> local_rows(numprocs, 0);
  local_rows[0] = 2;

  mat.reinit(MPI_COMM_WORLD, csp, local_rows, local_rows, myid);

  if (myid == 0)
    mat.add(0, 0, 1.0);


  mat.compress(VectorOperation::add);

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
