// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test
// LinearAlgebra::TpetraWrappers::SparseMatrix<double>::copy_from(other_matrix)
// for the case of a non-contiguous set of rows.

#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_sparsity_pattern.h>
#include <deal.II/lac/vector_operation.h>

#include <iostream>

#include "../tests.h"


void
test()
{
  const unsigned int MyPID   = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int NumProc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (MyPID == 0)
    deallog << "NumProc=" << NumProc << std::endl;

  // create non-contiguous index set for NumProc > 1
  dealii::IndexSet parallel_partitioning(NumProc * 2);

  // non-contiguous
  parallel_partitioning.add_index(MyPID);
  parallel_partitioning.add_index(NumProc + MyPID);

  // create sparsity pattern from parallel_partitioning

  // The sparsity pattern corresponds to a [FE_DGQ<1>(p=0)]^2 FESystem,
  // on a triangulation in which each MPI process owns 2 cells,
  // with reordered dofs by its components, such that the rows in the
  // final matrix are locally not in a contiguous set.

  dealii::LinearAlgebra::TpetraWrappers::SparsityPattern<MemorySpace::Default>
    sp_M(parallel_partitioning, MPI_COMM_WORLD, 2);

  sp_M.add(MyPID, MyPID);
  sp_M.add(MyPID, NumProc + MyPID);
  sp_M.add(NumProc + MyPID, MyPID);
  sp_M.add(NumProc + MyPID, NumProc + MyPID);

  sp_M.compress();

  // create matrix with dummy entries on the diagonal
  dealii::LinearAlgebra::TpetraWrappers::SparseMatrix<double,
                                                      MemorySpace::Default>
    M0;
  M0.reinit(sp_M);
  M0 = 0;

  for (const auto &i : parallel_partitioning)
    M0.set(i, i, dealii::numbers::PI);

  M0.compress(dealii::VectorOperation::insert);

  ////////////////////////////////////////////////////////////////////////
  // test ::add(TrilinosScalar, SparseMatrix)
  //

  dealii::LinearAlgebra::TpetraWrappers::SparseMatrix<double,
                                                      MemorySpace::Default>
    M1;
  M1.reinit(sp_M); // avoid deep copy
  M1 = 0;
  M1.copy_from(M0);

  // check
  for (const auto &i : parallel_partitioning)
    {
      const auto &el = M1.el(i, i);

      if (MyPID == 0)
        deallog << "i = " << i << " , j = " << i << " , el = " << el
                << std::endl;

      AssertThrow(el == dealii::numbers::PI, dealii::ExcInternalError());
    }

  if (MyPID == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  test();
}
