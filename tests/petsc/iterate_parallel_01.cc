// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Check that we can iterate over the elements of a parallel
// PETSc matrix. This was made more difficult than necessary
// because calling PETScWrappers::MatrixBase::end(row) may
// access a row that is not actually stored on the current
// processor, but is also not the one-past-the-end row that
// we internally tested.
//
// The test is lightly adapted from one provided by Feimi Yu

#include "../tests.h"
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/base/index_set.h>
#include <iostream>
#include <vector>

void
test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  deallog << numproc << std::endl;

  // set up disjoint index sets and some overlapping ghosted elements
  IndexSet local_active(numproc*2);
  local_active.add_range(myid*2,myid*2+2);
  IndexSet local_relevant= local_active;
  local_relevant.add_range(0,1);

  DynamicSparsityPattern csp (local_relevant);

  for (unsigned int i=0; i<2*numproc; ++i)
    if (local_relevant.is_element(i))
      csp.add(i,i);

  csp.add(0,1);


  PETScWrappers::MPI::SparseMatrix mat;
  mat.reinit (local_active, local_active, csp, MPI_COMM_WORLD);

  Assert(mat.n()==numproc*2, ExcInternalError());
  Assert(mat.m()==numproc*2, ExcInternalError());

  // set local values
  mat.set(myid*2,myid*2, 1.0+myid*2.0);
  mat.set(myid*2+1,myid*2+1, 2.0+myid*2.0);

  mat.compress(VectorOperation::insert);

  //////////////////////////////////////////////
  /////This is a test for the local matrix iterator
  //////////////////////////////////////////////
  unsigned int start_row = mat.local_range().first;
  unsigned int end_row = mat.local_range().second;
  for (auto r = start_row; r < end_row; ++r)
    {
      for (auto itr = mat.begin(r); itr != mat.end(r); ++itr)
        {
          deallog << itr->row() << ' ' << itr->column() << ' ' << itr->value()
                  << std::endl;
        }
    }
  /////////////////////////////////////////////

  if (myid==0)
    deallog << "OK" << std::endl;
}


int
main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  MPILogInitAll log;

  test ();
}
