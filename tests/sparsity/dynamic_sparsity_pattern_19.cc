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



// check SparsityTools::gather_sparsity_pattern()

#include <deal.II/lac/sparsity_tools.h>

#include "../tests.h"


void
test()
{
  const unsigned int N = 10;
  MPI_Comm           comm(MPI_COMM_WORLD);
  const unsigned int myid = Utilities::MPI::this_mpi_process(comm);

  IndexSet owned(N);
  if (myid == 0)
    owned.add_range(0, 3);
  else if (myid == 1)
    owned.add_range(3, 7);
  else
    owned.add_range(7, 10);

  IndexSet ghost_range(owned);
  if (myid == 0)
    {
      ghost_range.add_range(3, 5);
    }
  else if (myid == 1)
    {
      ghost_range.add_range(1, 3);
      ghost_range.add_range(7, 9);
    }
  else
    {
      ghost_range.add_range(4, 7);
    }

  DynamicSparsityPattern dsp;
  dsp.reinit(N, N);
  // hard-code some sparsity
  if (myid == 0)
    {
      for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 5 + i; ++j)
          dsp.add(i, j);

      // add to ghost to make sure they are reset
      dsp.add(3, 9);
    }
  else if (myid == 1)
    {
      for (unsigned int i = 3; i < 5; ++i)
        for (unsigned int j = 0; j < 5 + i; ++j)
          dsp.add(i, j);

      for (unsigned int i = 5; i < 7; ++i)
        for (unsigned int j = i - 4; j < 10; ++j)
          dsp.add(i, j);
    }
  else
    {
      for (unsigned int i = 7; i < 10; ++i)
        for (unsigned int j = i - 4; j < 10; ++j)
          dsp.add(i, j);
    }

  deallog << "Before gather_sparsity_pattern:" << std::endl;
  dsp.print(deallog.get_file_stream());

  SparsityTools::gather_sparsity_pattern(dsp, owned, comm, ghost_range);

  deallog << "After gather_sparsity_pattern:" << std::endl;
  dsp.print(deallog.get_file_stream());
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll log;

  test();
  return 0;
}
