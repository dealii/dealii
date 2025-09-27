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


// check MPI_InitializeFinalize::register_request

#include <deal.II/base/mpi.h>

#include "../tests.h"

void
unnecessary(MPI_Comm comm)
{
  // We register this request, but actually decide to wait at the end of this
  // function, which means that the request is guaranteed to be NULL
  // already. Make sure this works:
  static MPI_Request request = MPI_REQUEST_NULL;
  Utilities::MPI::MPI_InitFinalize::register_request(request);

  int ierr = MPI_Ibarrier(comm, &request);
  AssertThrowMPI(ierr);

  ierr = MPI_Wait(&request, MPI_STATUS_IGNORE);
  AssertThrowMPI(ierr);
}


void
test(MPI_Comm comm)
{
  // Here we execute an Ibarrier at the end and a wait at the
  // beginning. Without the barrier, multiple calls to this function fail as
  // messages are mixed up. With it, the algorithm works, but not registering
  // the request causes an error as well.

  static MPI_Request request = MPI_REQUEST_NULL;
  Utilities::MPI::MPI_InitFinalize::register_request(request);
  int ierr = MPI_Wait(&request, MPI_STATUS_IGNORE);
  AssertThrowMPI(ierr);

  int        tag     = 12345;
  const auto my_rank = Utilities::MPI::this_mpi_process(comm);
  const auto n_ranks = Utilities::MPI::n_mpi_processes(comm);

  if (my_rank == 0)
    {
      std::set<int> received_from;
      MPI_Status    status;

      for (unsigned int n = 1; n < n_ranks; ++n)
        {
          unsigned int value;
          MPI_Recv(&value, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, tag, comm, &status);

          AssertThrow(received_from.count(status.MPI_SOURCE) == 0,
                      ExcMessage("oh no!"));
          received_from.insert(status.MPI_SOURCE);
        }
    }
  else
    {
      unsigned int value = 123;
      int          dest  = 0;
      MPI_Send(&value, 1, MPI_UNSIGNED, dest, tag, comm);
    }

  ierr = MPI_Ibarrier(comm, &request);
  AssertThrowMPI(ierr);
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  mpi_initlog();

  unnecessary(MPI_COMM_WORLD);
  unnecessary(MPI_COMM_WORLD);

  test(MPI_COMM_WORLD);
  test(MPI_COMM_WORLD);
  test(MPI_COMM_WORLD);
}
