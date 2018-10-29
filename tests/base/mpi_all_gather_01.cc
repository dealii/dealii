// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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
  // Creating the local array of points
  std::vector<Point<3>> local_points(my_proc + 1);
  for (unsigned int i = 0; i < my_proc + 1; ++i)
    local_points[i] = Point<3>(my_proc, -my_proc, i);

  auto gathered_points =
    Utilities::MPI::all_gather(MPI_COMM_WORLD, local_points);
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
