// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
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
