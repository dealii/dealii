// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2017 by the deal.II authors
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
