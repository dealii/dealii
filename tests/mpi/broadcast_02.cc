// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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

// Test Utilities::MPI::broadcast() including sending a large message
// with >2^31 elements

#include <deal.II/base/mpi.h>

#include "../tests.h"

void
check1()
{
  const auto my_proc = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  int x[3];
  if (my_proc == 0)
    {
      x[0] = 1;
      x[1] = 2;
      x[2] = 3;
    }

  Utilities::MPI::broadcast(x, 3, 0, MPI_COMM_WORLD);
  AssertThrow(x[0] == 1, ExcInternalError());
  AssertThrow(x[1] == 2, ExcInternalError());
  AssertThrow(x[2] == 3, ExcInternalError());
  deallog << "OK" << std::endl;
}

void
check2(std::size_t count)
{
  const auto my_proc = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  std::vector<char> x(count, '?');
  if (my_proc == 0)
    {
      // a classical baby shout:
      std::fill(x.begin(), x.end(), 'A');
      x[count - 1] = '!';
    }

  Utilities::MPI::broadcast(x.data(), count, 0, MPI_COMM_WORLD);

  AssertThrow(x[0] == 'A', ExcInternalError());
  AssertThrow(x[count - 1] == '!', ExcInternalError());
  deallog << "OK" << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  check1();
  check2(1ULL << 3);
  check2(1ULL << 31);
  check2((1ULL << 32) + 5);
}
