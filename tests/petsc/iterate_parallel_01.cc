// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that we can iterate over the elements of a parallel
// PETSc matrix. This was made more difficult than necessary
// because calling PETScWrappers::MatrixBase::end(row) may
// access a row that is not actually stored on the current
// processor, but is also not the one-past-the-end row that
// we internally tested.
//
// The test is lightly adapted from one provided by Feimi Yu

#include <deal.II/base/index_set.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_sparse_matrix.h>

#include <iostream>
#include <vector>

#include "../tests.h"

void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  deallog << numproc << std::endl;

  // set up disjoint index sets and some overlapping ghosted elements
  IndexSet local_active(numproc * 2);
  local_active.add_range(myid * 2, myid * 2 + 2);
  IndexSet local_relevant = local_active;
  local_relevant.add_range(0, 1);

  DynamicSparsityPattern csp(local_relevant);

  for (unsigned int i = 0; i < 2 * numproc; ++i)
    if (local_relevant.is_element(i))
      csp.add(i, i);

  csp.add(0, 1);


  PETScWrappers::MPI::SparseMatrix mat;
  mat.reinit(local_active, local_active, csp, MPI_COMM_WORLD);

  Assert(mat.n() == numproc * 2, ExcInternalError());
  Assert(mat.m() == numproc * 2, ExcInternalError());

  // set local values
  mat.set(myid * 2, myid * 2, 1.0 + myid * 2.0);
  mat.set(myid * 2 + 1, myid * 2 + 1, 2.0 + myid * 2.0);

  mat.compress(VectorOperation::insert);

  //////////////////////////////////////////////
  /////This is a test for the local matrix iterator
  //////////////////////////////////////////////
  unsigned int start_row = mat.local_range().first;
  unsigned int end_row   = mat.local_range().second;
  for (auto r = start_row; r < end_row; ++r)
    {
      for (auto itr = mat.begin(r); itr != mat.end(r); ++itr)
        {
          deallog << itr->row() << ' ' << itr->column() << ' ' << itr->value()
                  << std::endl;
        }
    }
  /////////////////////////////////////////////

  if (myid == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test();
}
