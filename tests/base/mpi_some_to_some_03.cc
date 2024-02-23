// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test Utilities::MPI::some_to_some when sending to myself as well.

#include <deal.II/base/mpi.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/point.h>

#include <tuple>
#include <vector>

#include "../tests.h"


void
test()
{
  const auto n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const auto my_proc = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  std::map<unsigned int, std::vector<unsigned int>> to_send;

  // send to myself and to my next one
  to_send[my_proc].emplace_back(my_proc + 1);
  to_send[(my_proc + 1) % n_procs].emplace_back((my_proc + 1) * 10);

  const auto received = Utilities::MPI::some_to_some(MPI_COMM_WORLD, to_send);

  const auto original = Utilities::MPI::some_to_some(MPI_COMM_WORLD, received);

  deallog << "Sent                      : "
          << Patterns::Tools::to_string(to_send) << std::endl;
  deallog << "Received                  : "
          << Patterns::Tools::to_string(received) << std::endl;
  deallog << "Received(Received) == Sent: "
          << Patterns::Tools::to_string(original) << std::endl;

  // now check that to_send and original are the same
  if (original == to_send)
    deallog << "OK" << std::endl;
  else
    deallog << "Not OK" << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test();
}
