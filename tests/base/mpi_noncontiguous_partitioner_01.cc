// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


// Test Utilities::MPI::NoncontiguousPartitioner for non-contiguous index space.

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.h>

#include "../tests.h"


void
test(const MPI_Comm comm)
{
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

  Utilities::MPI::NoncontiguousPartitioner<double> vector(index_set_has,
                                                          index_set_want,
                                                          comm);

  AlignedVector<double> src(index_set_has.n_elements());
  AlignedVector<double> dst(index_set_want.n_elements());

  src[0] = Utilities::MPI::this_mpi_process(comm) * 100 + 1;

  vector.update_values(dst, src);

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

  {
    deallog.push("all");
    test(comm);
    deallog.pop();
  }
}
