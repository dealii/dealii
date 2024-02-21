// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test for the internal::send_and_receive function in the case of
// collective communication

#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>

#include <vector>

#include "../tests.h"

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  auto n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  auto my_proc = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  // Creating the map array of points to be sent
  std::map<unsigned int, std::vector<Point<2>>> m;
  for (unsigned int i = 0; i < n_procs; ++i)
    if (i != my_proc)
      m[i].push_back(Point<2>(my_proc, -2.0 * my_proc));

  auto received_pts = Utilities::MPI::some_to_some(MPI_COMM_WORLD, m);

  bool test_passed = true;
  for (const auto &pt : received_pts)
    if (std::abs(pt.first - pt.second[0][0]) > 1e-12 ||
        std::abs(2.0 * pt.first + pt.second[0][1]) > 1e-12)
      {
        test_passed = false;
        deallog << "Error with point " << pt.second[0] << " received from rank "
                << pt.first << std::endl;
      }
  if (test_passed)
    deallog << "Test: ok" << std::endl;
  else
    deallog << "Test: FAILED" << std::endl;
}
