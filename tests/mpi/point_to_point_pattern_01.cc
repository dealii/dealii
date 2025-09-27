// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check if point-to-point communication
// works correctly with MPI_Reduce_scatter_block
// implementation.
// Also check that the simplified point-to-point communication
// computation agrees with the original function.
// Namely, the size of
// compute_point_to_point_communication_pattern()
// returned vector should agree with the number returned
// by compute_point_to_point_communications()

#include <deal.II/base/utilities.h>

#include <algorithm>

#include "../tests.h"


void
test_mpi()
{
  Assert(Utilities::MPI::job_supports_mpi(), ExcInternalError());

  unsigned int       myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  // select a few destinations
  std::vector<unsigned int> destinations;
  for (unsigned int i = 0; i < 3 + myid / 3; ++i)
    if ((myid + 17 * i) % numprocs != myid)
      destinations.push_back((myid + 17 * i) % numprocs);

  if (myid == 0)
    {
      deallog << "Processor 0 wants to send to ";
      for (unsigned int i = 0; i < destinations.size(); ++i)
        deallog << destinations[i] << ' ';
      deallog << std::endl;

      for (unsigned int p = 1; p < numprocs; ++p)
        {
          MPI_Status   status;
          unsigned int size = 0;
          MPI_Recv(&size, 1, MPI_UNSIGNED, p, 0, MPI_COMM_WORLD, &status);

          std::vector<unsigned int> dest(size);
          MPI_Recv(&dest[0], size, MPI_UNSIGNED, p, 0, MPI_COMM_WORLD, &status);

          deallog << "Processor " << p << " wants to send to ";
          for (unsigned int i = 0; i < size; ++i)
            deallog << dest[i] << ' ';
          deallog << std::endl;
        }
    }
  else
    {
      unsigned int size = destinations.size();
      MPI_Send(&size, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&destinations[0], size, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
    }


  if (myid == 0)
    deallog << "Exchanging data..." << std::endl;

  std::vector<unsigned int> origins =
    Utilities::MPI::compute_point_to_point_communication_pattern(MPI_COMM_WORLD,
                                                                 destinations);

  const unsigned int n_origins =
    Utilities::MPI::compute_n_point_to_point_communications(MPI_COMM_WORLD,
                                                            destinations);

  if (origins.size() != n_origins)
    deallog << "Size mismatch!" << std::endl;

  std::sort(origins.begin(), origins.end());

  if (myid == 0)
    {
      deallog << "Processor 0 will receive from ";
      for (unsigned int i = 0; i < origins.size(); ++i)
        deallog << origins[i] << ' ';
      deallog << std::endl;

      for (unsigned int p = 1; p < numprocs; ++p)
        {
          MPI_Status   status;
          unsigned int size = 0;
          MPI_Recv(&size, 1, MPI_UNSIGNED, p, 0, MPI_COMM_WORLD, &status);

          std::vector<unsigned int> orig(size);
          MPI_Recv(&orig[0], size, MPI_UNSIGNED, p, 0, MPI_COMM_WORLD, &status);

          deallog << "Processor " << p << " will receive from ";
          for (unsigned int i = 0; i < size; ++i)
            deallog << orig[i] << ' ';
          deallog << std::endl;
        }
    }
  else
    {
      unsigned int size = origins.size();
      MPI_Send(&size, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&origins[0], size, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();

      deallog.push("mpi");
      test_mpi();
      deallog.pop();
    }
  else
    test_mpi();
}
