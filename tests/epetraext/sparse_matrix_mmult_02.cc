// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check TrilinosWrappers::SparseMatrix::mmult in parallel

#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <algorithm>
#include <random>

#include "../tests.h"


void
value_rank_0(TrilinosWrappers::SparseMatrix &a)
{
  a.set(0, 0, 0);
  a.set(0, 2, 2);
  a.set(0, 15, 15);
  a.set(1, 4, 5);
  a.set(1, 9, 10);
  a.set(1, 10, 11);
  a.set(2, 0, 2);
  a.set(2, 13, 30);
  a.set(3, 7, 10);
  a.set(3, 10, 13);
  a.set(3, 18, 21);
  a.set(4, 0, 4);
  a.set(4, 1, 5);
  a.set(4, 16, 20);
  a.set(5, 0, 5);
  a.set(5, 10, 15);
  a.set(5, 13, 18);
  a.set(6, 1, 7);
  a.set(6, 7, 13);
  a.set(6, 8, 14);
  a.set(7, 13, 20);
  a.set(7, 18, 25);
  a.set(7, 11, 18);
  a.set(8, 1, 9);
  a.set(8, 10, 18);
  a.set(8, 16, 24);
  a.set(9, 8, 17);
  a.set(9, 13, 22);
  a.set(9, 14, 23);
}

void
value_rank_1(TrilinosWrappers::SparseMatrix &a)
{
  a.set(10, 15, 25);
  a.set(10, 0, 10);
  a.set(10, 2, 12);
  a.set(11, 10, 21);
  a.set(11, 9, 20);
  a.set(11, 4, 15);
  a.set(12, 13, 50);
  a.set(12, 0, 12);
  a.set(13, 10, 23);
  a.set(13, 18, 31);
  a.set(13, 7, 20);
  a.set(14, 16, 30);
  a.set(14, 0, 14);
  a.set(14, 1, 15);
  a.set(15, 10, 25);
  a.set(15, 13, 28);
  a.set(15, 0, 15);
  a.set(16, 7, 23);
  a.set(16, 1, 17);
  a.set(16, 8, 24);
  a.set(17, 11, 28);
  a.set(17, 13, 30);
  a.set(17, 18, 35);
  a.set(18, 10, 28);
  a.set(18, 16, 34);
  a.set(18, 1, 19);
  a.set(19, 13, 32);
  a.set(19, 14, 33);
  a.set(19, 8, 27);
}

void
test()
{
  const unsigned int my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  const unsigned int n_local_rows      = 10;
  const unsigned int n_entries_per_row = 3;
  const unsigned int size              = n_local_rows * n_procs;

  // Compare the result of a*b in parallel with a*b in serial
  IndexSet           partitioning(size);
  const unsigned int local_offset(n_local_rows * my_rank);
  partitioning.add_range(local_offset, local_offset + n_local_rows);
  partitioning.compress();
  TrilinosWrappers::SparseMatrix a(partitioning);
  TrilinosWrappers::SparseMatrix b(partitioning);
  if (my_rank == 0)
    {
      value_rank_0(a);
    }
  if (my_rank == 1)
    {
      value_rank_1(a);
    }
  a.compress(VectorOperation::insert);
  b.copy_from(a);
  TrilinosWrappers::SparseMatrix c;
  a.mmult(c, b);


  TrilinosWrappers::SparseMatrix serial_a(size, size, n_entries_per_row);
  TrilinosWrappers::SparseMatrix serial_b(size, size, n_entries_per_row);
  value_rank_0(serial_a);
  value_rank_1(serial_a);
  serial_a.compress(VectorOperation::insert);
  serial_b.copy_from(serial_a);
  TrilinosWrappers::SparseMatrix serial_c;
  serial_a.mmult(serial_c, serial_b);

  AssertThrow(serial_c.n_nonzero_elements() == c.n_nonzero_elements(),
              ExcInternalError());
  AssertThrow(serial_c.l1_norm() == c.l1_norm(), ExcInternalError());
  AssertThrow(serial_c.linfty_norm() == c.linfty_norm(), ExcInternalError());
  AssertThrow(serial_c.frobenius_norm() == c.frobenius_norm(),
              ExcInternalError());

  if (my_rank == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  initlog();

  test();
}
