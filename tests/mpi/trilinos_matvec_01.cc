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



// Test whether TrilinosWrappers::SparseMatrix::vmult gives same result with
// Trilinos vector and parallel distributed vector

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

  const unsigned int n_rows = 3;
  const unsigned int n_cols = 4;

  IndexSet row_partitioning(n_rows);
  IndexSet col_partitioning(n_cols);

  if (n_procs == 1)
    {
      row_partitioning.add_range(0, n_rows);
      col_partitioning.add_range(0, n_cols);
    }
  else if (n_procs == 2)
    {
      // row_partitioning should be { [0, 2), [2, n_rows) }
      // col_partitioning should be { [0, 2), [2, n_cols) }
      // col_relevant_set should be { [0, 3), [1, n_cols) }
      if (my_id == 0)
        {
          row_partitioning.add_range(0, 2);
          col_partitioning.add_range(0, 2);
        }
      else if (my_id == 1)
        {
          row_partitioning.add_range(2, n_rows);
          col_partitioning.add_range(2, n_cols);
        }
    }
  else
    DEAL_II_NOT_IMPLEMENTED();

  TrilinosWrappers::SparsityPattern sp(row_partitioning,
                                       col_partitioning,
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
  A.clear();
  A.reinit(sp);
  if (my_id == 0)
    {
      A.add(0, 0, 1);
      A.add(0, 2, 1);
    }
  if ((n_procs == 1) || (my_id == 1))
    A.add(2, 3, 2.0);
  A.compress(VectorOperation::add);

  TrilinosWrappers::MPI::Vector x, y;
  x.reinit(col_partitioning, MPI_COMM_WORLD);
  y.reinit(row_partitioning, MPI_COMM_WORLD);

  LinearAlgebra::distributed::Vector<double> dx(col_partitioning,
                                                col_partitioning,
                                                MPI_COMM_WORLD),
    dy(row_partitioning, row_partitioning, MPI_COMM_WORLD);

  for (unsigned int i = 0; i < col_partitioning.n_elements(); ++i)
    {
      const unsigned int global_index = col_partitioning.nth_index_in_set(i);
      dx(global_index)                = random_value<double>();
      x(global_index)                 = dx(global_index);
    }
  dy = 1.;

  A.vmult(y, x);
  A.vmult(dy, dx);

  // compare whether we got the same result
  // (should be no roundoff difference)
  for (unsigned int i = 0; i < row_partitioning.n_elements(); ++i)
    {
      const unsigned int global_index = row_partitioning.nth_index_in_set(i);
      AssertThrow(dy(global_index) == y(global_index), ExcInternalError());
    }
  if (my_id == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  unsigned int       myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}
