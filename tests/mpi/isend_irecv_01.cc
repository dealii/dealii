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


// Test Utilities::MPI::isend() and Utilities::MPI::irecv().

#include <deal.II/base/mpi.h>

#include "../tests.h"



void
test(const MPI_Comm comm)
{
  using MessageType = unsigned int;

  // First post a bunch of receives before we start working on the
  // sends. We expect ten sends from each other process.
  std::multimap<unsigned int, Utilities::MPI::Future<MessageType>>
    receive_futures;
  for (unsigned int rank = 0; rank < Utilities::MPI::n_mpi_processes(comm);
       ++rank)
    if (rank != Utilities::MPI::this_mpi_process(comm))
      for (unsigned int message = 0; message < 10; ++message)
        {
          receive_futures.emplace(rank,
                                  Utilities::MPI::irecv<MessageType>(comm,
                                                                     rank));
        }

  // Wait a bit -- make sure all of those std::future-like objects are
  // actually delayed a bit before we query them below.
  sleep(10);

  // Send ten packages each to each of the other processes. Each
  // packet is simply an unsigned integer with that other process's
  // rank.
  std::vector<Utilities::MPI::Future<void>> send_futures;
  for (unsigned int rank = 0; rank < Utilities::MPI::n_mpi_processes(comm);
       ++rank)
    if (rank != Utilities::MPI::this_mpi_process(comm))
      for (unsigned int message = 0; message < 10; ++message)
        {
          send_futures.emplace_back(
            Utilities::MPI::isend(MessageType(rank), comm, rank));
        }

  // Finally wait for the receives to be done as well.
  for (auto &receive_future : receive_futures)
    {
      // Wait for the receive, obtain the message, and make sure that
      // it is what we expected:
      const MessageType message = receive_future.second.get();
      Assert(message == Utilities::MPI::this_mpi_process(comm),
             ExcMessage("On process " +
                        std::to_string(Utilities::MPI::this_mpi_process(comm)) +
                        ", we were expected a message from process " +
                        std::to_string(receive_future.first) +
                        " with a message value equal to the destination "
                        "process (that is, our own rank). But we got " +
                        std::to_string(message)));
      deallog << "Receive future returned successfully." << std::endl;
    }

  // Finally wait for the sends to be done as well.
  for (auto &send_future : send_futures)
    {
      send_future.wait();
      deallog << "Send future returned successfully." << std::endl;
    }
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const MPI_Comm comm = MPI_COMM_WORLD;

  test(comm);
}
