// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Trilinos PreconditionBlock* used to fail if one processor has no locally
// owned rows

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include "../tests.h"


template <class Prec>
void
test()
{
  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  // only proc 0 and 1 own rows
  IndexSet rows(4 * n_procs);
  if (myid == 0)
    rows.add_range(0, 2 * n_procs);
  else if (myid == 1)
    rows.add_range(2 * n_procs, 4 * n_procs);
  rows.compress();

  TrilinosWrappers::MPI::Vector src(rows, MPI_COMM_WORLD),
    dst(rows, MPI_COMM_WORLD);

  TrilinosWrappers::SparseMatrix mat(rows, rows, MPI_COMM_WORLD);
  for (const auto &row : rows)
    {
      const unsigned int i = row;
      mat.set(i, i, 100.);
      for (unsigned int j = 0; j < mat.n(); ++j)
        if (i != j)
          mat.set(i, j, i * j * .5 + .5);
    }

  for (unsigned int i = 0; i < src.size(); ++i)
    src(i) = i;

  mat.compress(VectorOperation::insert);
  src.compress(VectorOperation::insert);

  Prec                          preconditioner;
  typename Prec::AdditionalData data;
  data.block_size = 4;

  preconditioner.initialize(mat, data);
  preconditioner.vmult(dst, src);

  deallog << "dst: " << dst.l2_norm() << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<TrilinosWrappers::PreconditionBlockJacobi>();
  test<TrilinosWrappers::PreconditionBlockSOR>();
  test<TrilinosWrappers::PreconditionBlockSSOR>();
}
