// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
