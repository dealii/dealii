// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


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
