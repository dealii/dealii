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



// check DynamicSparsityPattern::gather()

#include <deal.II/lac/dynamic_sparsity_pattern.h>

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

  DynamicSparsityPattern dsp;

  dsp.reinit(N, N, owned, comm);
  // hard-code some sparsity
  if (myid == 0)
    {
      for (unsigned int i = 0; i <= 4; ++i)
        for (unsigned int j = 0; j <= 4; ++j)
          dsp.add(i, j);
    }
  else if (myid == 1)
    {
      for (unsigned int i = 1; i <= 4; ++i)
        for (unsigned int j = 1; j <= 4 + i; ++j)
          dsp.add(i, j);

      for (unsigned int i = 5; i <= 8; ++i)
        for (unsigned int j = i - 4; j <= 8; ++j)
          dsp.add(i, j);
    }
  else
    {
      for (unsigned int j = 4; j <= 8; ++j)
        {
          dsp.add(4, j);
          dsp.add(9, j + 1);
        }

      for (unsigned int i = 5; i <= 8; ++i)
        for (unsigned int j = 4; j <= 9; ++j)
          dsp.add(i, j);
    }

  deallog << "Before compress:" << std::endl;
  dsp.print(deallog.get_file_stream());
  dsp.compress();
  deallog << "After compress:" << std::endl;
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
