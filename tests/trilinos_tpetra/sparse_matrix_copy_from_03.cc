// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Like sparse_matrix_copy_from_02.cc, but with matrix columns that reference
// entries on other MPI processes. Also check matrix structure.

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
    sp_M(parallel_partitioning, MPI_COMM_WORLD, 3);

  sp_M.add(MyPID, MyPID);
  sp_M.add(MyPID, NumProc + MyPID);
  sp_M.add(NumProc + MyPID, MyPID);
  sp_M.add(NumProc + MyPID, NumProc + MyPID);

  // Add an entry that is not in the locally owned rows, to test that it is
  // properly included This corresponds to a contributions from a ghost cell to
  // the locally owned DoF
  sp_M.add(MyPID, (MyPID + 1) % NumProc);

  sp_M.compress();

  // create matrix with dummy entries on the diagonal
  dealii::LinearAlgebra::TpetraWrappers::SparseMatrix<double,
                                                      MemorySpace::Default>
    M0;
  M0.reinit(sp_M);
  M0 = 0;

  M0.set(MyPID, MyPID, dealii::numbers::PI);
  M0.set(MyPID, NumProc + MyPID, 1.0);
  M0.set(NumProc + MyPID, MyPID, 2.0);
  M0.set(NumProc + MyPID, NumProc + MyPID, dealii::numbers::PI);

  M0.set(MyPID, (MyPID + 1) % NumProc, 3.0);

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

  // check entries
  for (const auto &i : parallel_partitioning)
    {
      const auto &el = M1.el(i, i);

      if (MyPID == 0)
        deallog << "i = " << i << " , j = " << i << " , el = " << el
                << std::endl;

      AssertThrow(el == dealii::numbers::PI, dealii::ExcInternalError());
    }

  const auto &el = M1.el(MyPID, NumProc + MyPID);
  deallog << "i = " << MyPID << " , j = " << NumProc + MyPID << " , el = " << el
          << std::endl;
  AssertThrow(el == 1.0, dealii::ExcInternalError());

  const auto &el1 = M1.el(NumProc + MyPID, MyPID);
  deallog << "i = " << NumProc + MyPID << " , j = " << MyPID
          << " , el = " << el1 << std::endl;
  AssertThrow(el1 == 2.0, dealii::ExcInternalError());

  const auto &el2 = M1.el(MyPID, (MyPID + 1) % NumProc);
  deallog << "i = " << MyPID << " , j = " << (MyPID + 1) % NumProc
          << " , el = " << el2 << std::endl;
  AssertThrow(el2 == 3.0, dealii::ExcInternalError());

  // check other matrix properties, structure, etc.
  AssertThrow(M1.is_compressed() == M0.is_compressed(),
              dealii::ExcInternalError());
  AssertThrow(M1.m() == M0.m(), dealii::ExcInternalError());
  AssertThrow(M1.n() == M0.n(), dealii::ExcInternalError());
  AssertThrow(M1.local_size() == M0.local_size(), dealii::ExcInternalError());
  AssertThrow(M1.local_range() == M0.local_range(), dealii::ExcInternalError());
  AssertThrow(M1.n_nonzero_elements() == M0.n_nonzero_elements(),
              dealii::ExcInternalError());

  if (MyPID == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  mpi_initlog();
  test();
}
