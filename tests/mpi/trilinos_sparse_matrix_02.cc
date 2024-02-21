// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// This demonstrates a bug where copying a Trilinos matrix and then modifying
// the source will edit both matrices.

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const unsigned int my_id   = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  const unsigned int n_rows = 4;

  IndexSet locally_owned(n_rows);

  if (n_procs == 1)
    {
      locally_owned.add_range(0, n_rows);
    }
  else if (n_procs == 2)
    {
      // should be { [0, 2), [2, n_rows) }
      if (my_id == 0)
        locally_owned.add_range(0, 2);
      else if (my_id == 1)
        locally_owned.add_range(2, n_rows);
    }
  else
    DEAL_II_NOT_IMPLEMENTED();

  TrilinosWrappers::SparsityPattern sp(locally_owned,
                                       locally_owned,
                                       MPI_COMM_WORLD);
  if (my_id == 0)
    {
      sp.add(0, 0);
      sp.add(0, 2);
    }
  if ((n_procs == 1) || (my_id == 1))
    sp.add(2, 3);
  sp.compress();

  TrilinosWrappers::SparseMatrix A;
  A.reinit(sp);
  TrilinosWrappers::SparseMatrix B;
  B.reinit(A);

  A.add(0, 0, 0.1);
  A.add(0, 2, 0.2);
  if ((n_procs == 1) || (my_id == 1))
    A.add(2, 3, 0.3);

  double l1a = (n_procs == 1) ? 0.3 : 0.4;
  double l1b = n_procs * 1.2;

  A.compress(VectorOperation::add);
  deallog << "1: " << A.l1_norm() << ' ' << B.l1_norm() << " (should be " << l1a
          << " 0.0)" << std::endl;

  deallog << "set B=A..." << std::endl;

  B.copy_from(A);

  deallog << "2: " << A.l1_norm() << ' ' << B.l1_norm() << " (should be " << l1a
          << ' ' << l1a << ')' << std::endl;

  if (my_id == 0)
    {
      deallog << "A(0,0)=" << A(0, 0) << std::endl;
      deallog << "B(0,0)=" << B(0, 0) << std::endl;
    }

  deallog << "reassemble A..." << std::endl;

  A = 0;
  A.add(0, 0, -1.2);
  A.compress(VectorOperation::add);
  deallog << "3: " << A.l1_norm() << ' ' << B.l1_norm() << " (should be " << l1b
          << ' ' << l1a << ')' << std::endl;

  if (my_id == 0)
    {
      deallog << "A(0,0)=" << A(0, 0) << std::endl;
      deallog << "B(0,0)=" << B(0, 0) << std::endl;
    }

  if (my_id == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  test();
}
