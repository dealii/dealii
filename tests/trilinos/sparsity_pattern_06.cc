// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test accessing the entries of TrilinosWrappers::SparsityPattern
// after compress has been called.
// The sparsity pattern used is the same as in sparse_matrix_add_03.

#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  const unsigned int MyPID   = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int NumProc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (MyPID == 0)
    deallog << "NumProc=" << NumProc << std::endl;

  // create non-contiguous index set for NumProc > 3
  dealii::IndexSet parallel_partitioning(NumProc * 2);

  // non-contiguous
  parallel_partitioning.add_index(MyPID);
  parallel_partitioning.add_index(NumProc + MyPID);

  // create sparsity pattern from parallel_partitioning

  // The sparsity pattern corresponds to a [FE_DGQ<1>(p=0)]^2 FESystem,
  // on a periodic triangulation in which each MPI process owns 2 cells,
  // with reordered dofs by its components, such that the rows in the
  // final matrix are locally not in a contiguous set.

  dealii::TrilinosWrappers::SparsityPattern sp(parallel_partitioning,
                                               MPI_COMM_WORLD,
                                               2);

  sp.add(MyPID, (NumProc + MyPID - 1) % NumProc);
  sp.add(MyPID, MyPID);
  sp.add(MyPID, (MyPID + 1) % NumProc);
  sp.add(MyPID + NumProc, (NumProc + MyPID - 1) % NumProc + NumProc);
  sp.add(MyPID + NumProc, MyPID + NumProc);
  sp.add(MyPID + NumProc, (MyPID + 1) % NumProc + NumProc);

  deallog << "before compress:" << std::endl;

  for (const auto &el : sp)
    {
      deallog << "index: " << el.index() << " = " << '(' << el.row() << " , "
              << el.column() << ')' << std::endl;
    }

  sp.compress();

  deallog << "after compress:" << std::endl;
  for (const auto &el : sp)
    {
      deallog << "index: " << el.index() << " = " << '(' << el.row() << " , "
              << el.column() << ')' << std::endl;
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  MPILogInitAll mpi_log;

  test();
}
