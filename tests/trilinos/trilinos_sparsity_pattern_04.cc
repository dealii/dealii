// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test several reinitialization options with a Trilinos sparsity pattern and
// Trilinos sparse matrix

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  TrilinosWrappers::SparsityPattern sp;

  deallog << "Creating entries..." << std::endl;

  IndexSet rows(2 * n_procs);
  rows.add_range(2 * myid, 2 * myid + 2);
  rows.compress();
  IndexSet writable_rows(rows);
  writable_rows.add_index(1);
  IndexSet columns(3 * n_procs);
  columns.add_range(3 * myid, 3 * myid + 3);
  columns.compress();

  // this creates a matrix with optimized path for off-processor entries
  sp.reinit(rows, columns, writable_rows, MPI_COMM_WORLD, 0u);

  for (unsigned int i = 2 * myid; i < 2 * myid + 2; ++i)
    for (unsigned int j = 0; j < 3 * n_procs; ++j)
      if ((i + 2 * j + 1) % 3 == 0)
        sp.add(i, j);

  sp.add(1, 0);
  sp.add(1, 4);

  sp.compress();

  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;
  deallog << "Number of entries: " << sp.n_nonzero_elements() << std::endl;
  deallog << "Number of rows: " << sp.n_rows() << std::endl;
  deallog << "Number of columns: " << sp.n_cols() << std::endl;

  {
    // create matrix from sparsity pattern, add a few entries and output the
    // matrix norm
    TrilinosWrappers::SparseMatrix matrix;
    matrix.reinit(sp);
    double c = 0;
    for (unsigned int i = 2 * myid; i < 2 * myid + 2; ++i)
      for (unsigned int j = 0; j < 3 * n_procs; ++j)
        if ((i + 2 * j + 1) % 3 == 0)
          matrix.add(i, j, c++);

    matrix.add(1, 0, c++);
    matrix.add(1, 4, c++);

    matrix.compress(VectorOperation::add);

    deallog << "Matrix norm: " << matrix.frobenius_norm() << std::endl;
  }

  // reinit, this time without giving writable rows -> Trilinos must manage
  // ghost entries and throw away the off-processor sparsity pattern
  sp.reinit(rows, columns, MPI_COMM_WORLD, 0u);

  for (unsigned int i = 2 * myid; i < 2 * myid + 2; ++i)
    for (unsigned int j = 0; j < 3 * n_procs; ++j)
      if ((i + 2 * j + 1) % 3 == 0)
        sp.add(i, j);

  sp.add(1, 0);
  sp.add(1, 4);
  sp.add(2, 1);
  sp.add(2, 5);

  sp.compress();
  {
    // create matrix from sparsity pattern, add a few entries and output the
    // matrix norm
    TrilinosWrappers::SparseMatrix matrix;
    matrix.reinit(sp);
    double c = 0;
    for (unsigned int i = 2 * myid; i < 2 * myid + 2; ++i)
      for (unsigned int j = 0; j < 3 * n_procs; ++j)
        if ((i + 2 * j + 1) % 3 == 0)
          matrix.add(i, j, c++);

    matrix.add(1, 0, c++);
    matrix.add(1, 4, c++);
    matrix.add(2, 1, c++);
    matrix.add(2, 5, c++);

    matrix.compress(VectorOperation::add);

    deallog << "Matrix norm: " << matrix.frobenius_norm() << std::endl;
  }


  // now create again a pattern with writable rows
  sp.reinit(rows, columns, writable_rows, MPI_COMM_WORLD, 0u);

  for (unsigned int i = 2 * myid; i < 2 * myid + 2; ++i)
    for (unsigned int j = 0; j < 3 * n_procs; ++j)
      if ((i + 2 * j + 1) % 3 == 0)
        sp.add(i, j);

  sp.add(1, 0);
  sp.add(1, 4);

  sp.compress();

  {
    // create matrix from sparsity pattern, add a few entries and output the
    // matrix norm
    TrilinosWrappers::SparseMatrix matrix;
    matrix.reinit(sp);
    double c = 0;
    for (unsigned int i = 2 * myid; i < 2 * myid + 2; ++i)
      for (unsigned int j = 0; j < 3 * n_procs; ++j)
        if ((i + 2 * j + 1) % 3 == 0)
          matrix.add(i, j, c++);

    matrix.add(1, 0, c++);
    matrix.add(1, 4, c++);

    matrix.compress(VectorOperation::add);

    deallog << "Matrix norm: " << matrix.frobenius_norm() << std::endl;
  }


  // finally, initialize the Trilinos sparsity pattern from a
  // DynamicSparsityPattern
  writable_rows.add_index(2);
  writable_rows.compress();
  DynamicSparsityPattern dsp(rows.size(), columns.size(), writable_rows);
  for (unsigned int i = 2 * myid; i < 2 * myid + 2; ++i)
    for (unsigned int j = 0; j < 3 * n_procs; ++j)
      if ((i + 2 * j + 1) % 3 == 0)
        dsp.add(i, j);

  dsp.add(1, 0);
  dsp.add(1, 4);
  dsp.add(2, 1);
  dsp.add(2, 5);

  sp.reinit(rows, columns, dsp, MPI_COMM_WORLD, true);
  sp.compress();
  {
    // create matrix from sparsity pattern, add a few entries and output the
    // matrix norm
    TrilinosWrappers::SparseMatrix matrix;
    matrix.reinit(sp);
    double c = 0;
    for (unsigned int i = 2 * myid; i < 2 * myid + 2; ++i)
      for (unsigned int j = 0; j < 3 * n_procs; ++j)
        if ((i + 2 * j + 1) % 3 == 0)
          matrix.add(i, j, c++);

    matrix.add(1, 0, c++);
    matrix.add(1, 4, c++);
    matrix.add(2, 1, c++);
    matrix.add(2, 5, c++);

    matrix.compress(VectorOperation::add);

    deallog << "Matrix norm: " << matrix.frobenius_norm() << std::endl;
  }

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();

      test();
    }
  else
    {
      test();
    }
}
