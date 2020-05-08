// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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


// Test Utilities::MPI::NoncontiguousPartitioner for padding.

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.h>

#include "../tests.h"


void
test(const MPI_Comm                       comm,
     std::vector<types::global_dof_index> index_set_has,
     std::vector<types::global_dof_index> index_set_want)
{
  Utilities::MPI::NoncontiguousPartitioner vector;
  vector.reinit(index_set_has, index_set_want, comm);

  AlignedVector<double> src(index_set_has.size(), 0);
  AlignedVector<double> dst(index_set_want.size(), 0);

  for (unsigned int i = 0; i < index_set_has.size(); i++)
    src[i] = Utilities::MPI::this_mpi_process(comm) * 100 + i;

  vector.export_to_ghosted_array(ArrayView<const double>(src.data(),
                                                         src.size()),
                                 ArrayView<double>(dst.data(), dst.size()));

  for (size_t i = 0; i < src.size(); i++)
    deallog << static_cast<int>(src[i]) << " ";
  deallog << std::endl;
  for (size_t i = 0; i < dst.size(); i++)
    deallog << static_cast<int>(dst[i]) << " ";
  deallog << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const MPI_Comm comm = MPI_COMM_WORLD;

  const unsigned int rank = Utilities::MPI::this_mpi_process(comm);

  {
    deallog.push("padding-non");

    if (rank == 0)
      test(comm, {0, 1, 2, 3}, {4, 5, 6, 7});
    else
      test(comm, {4, 5, 6, 7}, {0, 1, 2, 3});
    deallog.pop();
  }

  {
    deallog.push("padding-src");

    if (rank == 0)
      test(comm, {0, 1, numbers::invalid_dof_index, 2, 3}, {4, 5, 6, 7});
    else
      test(comm, {4, 5, 6, 7}, {0, 1, 2, 3});
    deallog.pop();
  }

  {
    deallog.push("padding-dst");

    if (rank == 0)
      test(comm, {0, 1, 2, 3}, {4, 5, numbers::invalid_dof_index, 6, 7});
    else
      test(comm, {4, 5, 6, 7}, {0, 1, 2, 3});
    deallog.pop();
  }
}
