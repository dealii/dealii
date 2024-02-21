// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test for Utilities::MPI::scatter

#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>

#include <vector>

#include "../tests.h"

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  const unsigned int root_process = 0;
  const auto         n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const auto         my_proc = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  // Creating the local array of points
  std::vector<Point<3>> local_points(my_proc + 1);
  for (unsigned int i = 0; i < my_proc + 1; ++i)
    local_points[i] = Point<3>(my_proc, -my_proc, i);

  // send to process 0
  const auto gathered_points =
    Utilities::MPI::gather(MPI_COMM_WORLD, local_points, root_process);

  // scatter from process 0
  const auto scattered_points =
    Utilities::MPI::scatter(MPI_COMM_WORLD, gathered_points, root_process);

  int test_passed = 1;

  if (scattered_points.size() != my_proc + 1)
    {
      test_passed = 0;
      deallog << "Error: Points on rank " << my_proc << " have wrong size. "
              << std::endl;
    }
  for (unsigned int p = 0; p < scattered_points.size(); ++p)
    if (scattered_points[p][0] != (double)my_proc ||
        scattered_points[p][1] != (double)-my_proc ||
        scattered_points[p][2] != (double)p)
      {
        test_passed = 0;
        deallog << "Error with point " << p << " on rank " << my_proc
                << std::endl;
      }

  if (Utilities::MPI::min(test_passed, MPI_COMM_WORLD))
    deallog << "Test: ok" << std::endl;
  else
    deallog << "Test: FAILED" << std::endl;
}
