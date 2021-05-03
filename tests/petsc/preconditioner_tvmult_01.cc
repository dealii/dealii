// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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


// Test Tvmult of PreconditionBase
// This test is based on reinit_preconditioner_01

#include <deal.II/base/index_set.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_sparse_matrix.h>

#include <iostream>
#include <vector>

#include "../tests.h"

template <class PRE>
void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  deallog << numproc << std::endl;

  // each processor owns 2 indices and all
  // are ghosting Element 1 (the second)

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

  {
    PETScWrappers::MPI::Vector src, dst;
    src.reinit(local_active, MPI_COMM_WORLD);
    dst.reinit(local_active, MPI_COMM_WORLD);
    src(myid * 2) = 1.0;
    src.compress(VectorOperation::insert);

    PRE pre;
    pre.initialize(mat);
    pre.Tvmult(dst, src);
    dst.print(deallog.get_file_stream());

    mat.add(0, 0, 1.0);
    mat.compress(VectorOperation::add);

    pre.initialize(mat);
    pre.Tvmult(dst, src);
    dst.print(deallog.get_file_stream());
  }

  if (myid == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<PETScWrappers::PreconditionJacobi>();
  test<PETScWrappers::PreconditionBlockJacobi>();
  test<PETScWrappers::PreconditionBoomerAMG>();
  test<PETScWrappers::PreconditionNone>();
}
