// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test for Utilities::MPI::gather

#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>

#include <vector>

#include "../tests.h"

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  unsigned int root_process = 0;
  auto         n_procs      = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  auto         my_proc      = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  // Creating the local array of points
  std::vector<Point<3>> local_points(my_proc + 1);
  for (unsigned int i = 0; i < my_proc + 1; ++i)
    local_points[i] = Point<3>(my_proc, -my_proc, i);

  auto gathered_points =
    Utilities::MPI::gather(MPI_COMM_WORLD, local_points, root_process);

  if (my_proc == root_process)
    {
      bool test_passed = true;
      for (unsigned int i = 0; i < n_procs; ++i)
        {
          if (gathered_points[i].size() != i + 1)
            {
              test_passed = false;
              deallog << "Error: Points received from rank " << i
                      << " have wrong size. " << std::endl;
            }
          for (unsigned int p = 0; p < gathered_points[i].size(); ++p)
            if (gathered_points[i][p][0] != (double)i ||
                gathered_points[i][p][1] != (double)-i ||
                gathered_points[i][p][2] != (double)p)
              {
                test_passed = false;
                deallog << "Error with point " << p << " from rank " << i
                        << std::endl;
              }
        }
      if (test_passed)
        deallog << "Test: ok" << std::endl;
      else
        deallog << "Test: FAILED" << std::endl;
    }
}
