// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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



// Tests IndexSet retrieval of Trilinos sparsity patterns

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
  IndexSet columns(3 * n_procs);
  columns.add_range(3 * myid, 3 * myid + 3);
  columns.compress();

  sp.reinit(rows, columns, MPI_COMM_WORLD, 0u);
  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;

  for (unsigned int i = 2 * myid; i < 2 * myid + 2; ++i)
    for (unsigned int j = 0; j < 3 * n_procs; ++j)
      if ((i + 2 * j + 1) % 3 == 0)
        sp.add(i, j);

  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;

  sp.compress();

  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;
  deallog << "Number of entries: " << sp.n_nonzero_elements() << std::endl;
  deallog << "Number of rows: " << sp.n_rows() << std::endl;
  deallog << "Number of columns: " << sp.n_cols() << std::endl;

  deallog << "Checks: ";
  IndexSet stored_rows = sp.locally_owned_range_indices();
  IndexSet stored_cols = sp.locally_owned_domain_indices();
  AssertThrow(stored_rows == rows, ExcInternalError());
  AssertThrow(stored_cols == columns, ExcInternalError());

  const unsigned int stored_n_procs =
    Utilities::MPI::n_mpi_processes(sp.get_mpi_communicator());
  const unsigned int stored_myid =
    Utilities::MPI::this_mpi_process(sp.get_mpi_communicator());

  AssertThrow(stored_n_procs == n_procs, ExcInternalError());
  AssertThrow(stored_myid == myid, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

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
