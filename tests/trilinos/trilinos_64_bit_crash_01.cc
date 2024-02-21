// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// extracted from deal_solver_02. at the time of writing this test, we ran
// into weird crashes with Trilinos when in 64 bit mode
//
// this has been fixed in Trilinos versions after 11.0.3. see
//   https://software.sandia.gov/bugzilla/show_bug.cgi?id=5884


#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include "../tests.h"


template <typename int_type>
void
test()
{
  const Epetra_Map map(int_type(1), 0, Epetra_MpiComm(MPI_COMM_SELF));

  int             n_entries_per_row[1] = {1};
  Epetra_CrsGraph graph(Copy, map, map, &n_entries_per_row[0], true);

  int_type row_indices[1] = {0};
  graph.InsertGlobalIndices(int_type(0), 1, &row_indices[0]);

  graph.FillComplete(map, map);
}


int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(4);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  deallog << "32bit" << std::endl;
  test<int>();

  deallog << "64bit" << std::endl;
  test<long long int>();

  deallog << "OK" << std::endl;
}
