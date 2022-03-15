// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2018 by the deal.II authors
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



// Test whether TrilinosWrappers::SparseMatrix::(T)vmult(_add) gives same
// result with Trilinos vector and parallel distributed vector when one
// processor does not have any rows

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
          row_partitioning.add_range(0, n_rows);
          col_partitioning.add_range(0, n_cols);
        }
    }
  else
    Assert(false, ExcNotImplemented());

  TrilinosWrappers::SparsityPattern sp(row_partitioning,
                                       col_partitioning,
                                       MPI_COMM_WORLD);
  if (my_id == 0)
    {
      sp.add(0, 0);
      sp.add(0, 2);
      sp.add(2, 3);
    }
  sp.compress();

  TrilinosWrappers::SparseMatrix A;
  A.clear();
  A.reinit(sp);
  if (my_id == 0)
    {
      A.add(0, 0, 1);
      A.add(0, 2, 1);
      A.add(2, 3, 2.0);
    }
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

  A.vmult_add(y, x);
  A.vmult_add(dy, dx);

  // compare whether we got the same result
  // (should be no roundoff difference)
  for (unsigned int i = 0; i < row_partitioning.n_elements(); ++i)
    {
      const unsigned int global_index = row_partitioning.nth_index_in_set(i);
      AssertThrow(dy(global_index) == y(global_index), ExcInternalError());
    }

  A.Tvmult(x, y);
  A.Tvmult(dx, dy);
  for (unsigned int i = 0; i < col_partitioning.n_elements(); ++i)
    {
      const unsigned int global_index = col_partitioning.nth_index_in_set(i);
      AssertThrow(dx(global_index) == x(global_index), ExcInternalError());
    }

  A.Tvmult_add(x, y);
  A.Tvmult_add(dx, dy);
  for (unsigned int i = 0; i < col_partitioning.n_elements(); ++i)
    {
      const unsigned int global_index = col_partitioning.nth_index_in_set(i);
      AssertThrow(dx(global_index) == x(global_index), ExcInternalError());
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
