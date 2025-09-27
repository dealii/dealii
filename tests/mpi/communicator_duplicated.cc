// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check MPI::duplicate_communicator, free_communicator, and
// DuplicatedCommunicator


#include <deal.II/base/mpi.h>

#include "../tests.h"


void
test(const MPI_Comm comm)
{
  int tag = 12345;

  const auto my_rank = Utilities::MPI::this_mpi_process(comm);

  MPI_Comm comm2 = Utilities::MPI::duplicate_communicator(comm);
  Utilities::MPI::DuplicatedCommunicator pcomm3(comm);

  if (my_rank == 1)
    {
      int         value[3] = {1, 2, 3};
      int         dest     = 0;
      MPI_Request requests[3];

      MPI_Isend(&value[0], 1, MPI_UNSIGNED, dest, tag, comm, &requests[0]);
      MPI_Isend(&value[1], 1, MPI_UNSIGNED, dest, tag, comm2, &requests[1]);
      MPI_Isend(&value[2], 1, MPI_UNSIGNED, dest, tag, *pcomm3, &requests[2]);
      MPI_Waitall(3, requests, MPI_STATUSES_IGNORE);
    }

  if (my_rank == 0)
    {
      int value[3];
      int src = 1;

      // receive in reverse order, if duplication of communicators worked,
      // these won't be mixed up!
      MPI_Recv(
        &value[2], 1, MPI_UNSIGNED, src, tag, *pcomm3, MPI_STATUS_IGNORE);
      MPI_Recv(&value[1], 1, MPI_UNSIGNED, src, tag, comm2, MPI_STATUS_IGNORE);
      MPI_Recv(&value[0], 1, MPI_UNSIGNED, src, tag, comm, MPI_STATUS_IGNORE);

      deallog << value[0] << ' ' << value[1] << ' ' << value[2] << std::endl;
    }

  Utilities::MPI::free_communicator(comm2);
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  mpi_initlog();

  test(MPI_COMM_WORLD);
}
