// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


// Test Utilities::MPI::NoncontiguousPartitioner for variable data transfer.

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.templates.h>

#include <vector>

#include "../tests.h"


void
test(const MPI_Comm comm)
{
  // prepare index sets
  IndexSet index_set_has(4);
  IndexSet index_set_want(4);

  if (Utilities::MPI::this_mpi_process(comm) == 0)
    {
      index_set_has.add_index(1);
      index_set_want.add_index(2);
    }
  else
    {
      index_set_has.add_index(2);
      index_set_want.add_index(1);
      index_set_want.add_index(2);
    }

  Utilities::MPI::NoncontiguousPartitioner vector(index_set_has,
                                                  index_set_want,
                                                  comm);

  // prepare containers and data
  AlignedVector<std::vector<int>> src(index_set_has.n_elements());
  AlignedVector<std::vector<int>> dst(index_set_want.n_elements());

  src[0] = std::vector<int>(Utilities::MPI::this_mpi_process(comm) + 1);

  // exchange data
  std::map<unsigned int, std::vector<char>> buffers;
  std::vector<unsigned int> sizes(vector.temporary_storage_size(), 0);

  vector.export_to_ghosted_array<std::vector<int>>(
    0,
    ArrayView<const std::vector<int>>(src.data(), src.size()),
    buffers,
    sizes,
    ArrayView<std::vector<int>>(dst.data(), dst.size()));

  // verify that sent data matches
  for (size_t i = 0; i < src.size(); i++)
    deallog << src[i].size() << " ";
  deallog << std::endl;
  for (size_t i = 0; i < dst.size(); i++)
    deallog << dst[i].size() << " ";
  deallog << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const MPI_Comm comm = MPI_COMM_WORLD;

  {
    deallog.push("all");
    test(comm);
    deallog.pop();
  }
}
