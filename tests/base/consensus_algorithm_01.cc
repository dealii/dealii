// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// Test ConsensusAlgorithms::AnonymousProcess.

#include <deal.II/base/mpi_consensus_algorithms.h>

#include "../tests.h"


void
test(const MPI_Comm &comm)
{
  const unsigned int my_rank = dealii::Utilities::MPI::this_mpi_process(comm);
  const unsigned int n_rank  = dealii::Utilities::MPI::n_mpi_processes(comm);

  using T1 = unsigned int;
  using T2 = unsigned int;

  dealii::Utilities::MPI::ConsensusAlgorithms::AnonymousProcess<T1, T2> process(
    [&]() {
      std::vector<unsigned int> result{(my_rank + 1) % n_rank};
      return result;
    },
    [&](const unsigned int other_rank, std::vector<T1> &send_buffer) {
      send_buffer.push_back(my_rank);
    },
    [&](const unsigned int &   other_rank,
        const std::vector<T1> &buffer_recv,
        std::vector<T2> &      request_buffer) {
      AssertDimension(other_rank, buffer_recv.front());
      deallog << "ConsensusAlgorithmProcess::answer_request() passed!"
              << std::endl;
      request_buffer.push_back(my_rank);
    },
    [&](const unsigned int other_rank, std::vector<T2> &recv_buffer) {
      recv_buffer.resize(1);
    },
    [&](const unsigned int other_rank, const std::vector<T2> &recv_buffer) {
      AssertDimension(other_rank, recv_buffer.front());
      deallog << "ConsensusAlgorithmProcess::function_read_answer() passed!"
              << std::endl;
    });
  dealii::Utilities::MPI::ConsensusAlgorithms::Selector<T1, T2>(process, comm)
    .run();
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const MPI_Comm comm = MPI_COMM_WORLD;

  test(comm);
}
