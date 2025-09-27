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



// test TrilinosWrappers::SparseMatrix::iterator semantics. make sure
// that rows not stored locally look like they're empty
//
// this test is like _04 but uses the const_iterator type which goes
// through a different code path in one place

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
  A.reinit(sp);
  if (my_id == 0)
    {
      A.set(0, 0, 0.1);
      A.set(0, 2, 0.2);
    }
  if ((n_procs == 1) || (my_id == 1))
    {
      A.set(0, 0, 0.1);
      A.set(0, 2, 0.2);
      A.set(2, 3, 0.3);
    }

  A.compress(VectorOperation::insert);

  const TrilinosWrappers::SparseMatrix &B = A;

  // now access elements by iterator. ensure that we can iterate over
  // all rows but that iterators into rows not stored locally just
  // look empty
  for (TrilinosWrappers::SparseMatrix::const_iterator p = B.begin();
       p != B.end();
       ++p)
    if (my_id == 0)
      {
        deallog << "Looking at entry (" << p->row() << ',' << p->column()
                << ") with value " << p->value() << std::endl;
        AssertThrow(p->row() == 0, ExcInternalError());
      }
    else
      {
        AssertThrow(p->row() == 2, ExcInternalError());
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
