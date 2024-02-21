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



// check setting off-processor entries of Epetra_CrsMatrix. It turns out that
// the underlying Epetra data structures actually add the entries even though
// we want to insert them only.

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <iostream>

#include "../tests.h"


void
test()
{
  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  IndexSet rows(5 * n_procs);
  rows.add_range(5 * myid, 5 * (myid + 1));
  rows.compress();
  IndexSet columns(5);
  columns.add_range(5 * myid / n_procs, 5 * (myid + 1) / n_procs);
  columns.compress();

  {
    TrilinosWrappers::SparseMatrix m(rows, columns, MPI_COMM_WORLD);
    for (unsigned int i = 5 * myid; i < 5 * (myid + 1); ++i)
      for (unsigned int j = 0; j < m.n(); ++j)
        if ((i + 2 * j + 1) % 3 == 0)
          m.set(i, j, i * j * .5 + .5);

    m.compress(VectorOperation::insert);
    deallog << "Matrix norm set local once: " << m.frobenius_norm()
            << std::endl;
  }

  {
    TrilinosWrappers::SparseMatrix m(rows, columns, MPI_COMM_WORLD);
    for (unsigned int i = 5 * myid; i < 5 * (myid + 1); ++i)
      for (unsigned int j = 0; j < m.n(); ++j)
        if ((i + 2 * j + 1) % 3 == 0)
          m.set(i, j, i * j * .5 + .5);

    for (unsigned int i = 5 * myid; i < 5 * (myid + 1); ++i)
      for (unsigned int j = 0; j < m.n(); ++j)
        if ((i + 2 * j + 1) % 3 == 0)
          m.set(i, j, i * j + 1.);

    m.compress(VectorOperation::insert);
    deallog << "Matrix norm set local twice: " << m.frobenius_norm()
            << std::endl;
  }

  {
    TrilinosWrappers::SparseMatrix m(rows, columns, MPI_COMM_WORLD);
    for (unsigned int i = 5 * myid; i < 5 * (myid + 1); ++i)
      for (unsigned int j = 0; j < m.n(); ++j)
        if ((i + 2 * j + 1) % 3 == 0)
          m.set(i, j, i * j * .5 + .5);

    if (myid == 1)
      m.set(1, 3, 10.);

    m.compress(VectorOperation::insert);
    deallog << "Matrix norm set non-local once: " << m.frobenius_norm()
            << std::endl;
  }

  {
    TrilinosWrappers::SparseMatrix m(rows, columns, MPI_COMM_WORLD);
    for (unsigned int i = 5 * myid; i < 5 * (myid + 1); ++i)
      for (unsigned int j = 0; j < m.n(); ++j)
        if ((i + 2 * j + 1) % 3 == 0)
          m.set(i, j, i * j * .5 + .5);

    m.set(1, 3, 10.);

    m.compress(VectorOperation::insert);
    deallog << "Matrix norm set non-local twice: " << m.frobenius_norm()
            << std::endl;
  }

  {
    TrilinosWrappers::SparseMatrix m(rows, columns, MPI_COMM_WORLD);
    for (unsigned int i = 5 * myid; i < 5 * (myid + 1); ++i)
      for (unsigned int j = 0; j < m.n(); ++j)
        if ((i + 2 * j + 1) % 3 == 0)
          m.set(i, j, i * j * .5 + .5);

    m.set(1, 3, 10.);
    m.set(2, 3, 2 * 3 * 0.5 + 0.5);

    m.compress(VectorOperation::insert);
    deallog << "Matrix norm set non-local twice: " << m.frobenius_norm()
            << std::endl;

    for (unsigned int i = 5 * myid; i < 5 * (myid + 1); ++i)
      for (unsigned int j = 0; j < m.n(); ++j)
        if ((i + 2 * j + 1) % 3 == 0)
          m.set(i, j, i * j * .5 + .5);

    m.set(1, 3, 10.);
    m.set(2, 3, 2 * 3 * 0.5 + 0.5);

    m.compress(VectorOperation::insert);
    deallog << "Matrix norm set twice, 2nd attempt: " << m.frobenius_norm()
            << std::endl;

    for (unsigned int i = 5 * myid; i < 5 * (myid + 1); ++i)
      for (unsigned int j = 0; j < m.n(); ++j)
        if ((i + 2 * j + 1) % 3 == 0)
          m.set(i, j, i * j * .5 + .5);

    m.set(1, 3, 10.);
    m.set(2, 3, 2 * 3 * 0.5 + 0.5);
    m.set(2, 3, 2 * 3 * 0.5 + 0.5);

    m.compress(VectorOperation::insert);
    deallog << "Matrix norm set twice-twice: " << m.frobenius_norm()
            << std::endl;
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
