// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
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


// Test ConsensusAlgorithms::selection().

#include <deal.II/base/mpi_consensus_algorithms.h>

#include "../tests.h"


void
test(const MPI_Comm &comm)
{
  const unsigned int my_rank = dealii::Utilities::MPI::this_mpi_process(comm);
  const unsigned int n_rank  = dealii::Utilities::MPI::n_mpi_processes(comm);

  using T1 = std::vector<unsigned int>;
  using T2 = std::vector<unsigned int>;

  const auto sources =
    dealii::Utilities::MPI::ConsensusAlgorithms::selector<T1, T2>(
      /* target_processes: */
      std::vector<unsigned int>{(my_rank + 1) % n_rank},
      /* create_request: */
      [my_rank](const unsigned int) { return T1({my_rank}); },
      /* answer_request: */
      [my_rank](const unsigned int other_rank, const T1 &request) {
        AssertDimension(other_rank, request.front());
        deallog << "ConsensusAlgorithmProcess::answer_request() passed!"
                << std::endl;
        return T2({my_rank});
      },
      /* process_answer: */
      [](const unsigned int other_rank, const T2 &answer) {
        AssertDimension(other_rank, answer.front());
        deallog << "ConsensusAlgorithmProcess::function_read_answer() passed!"
                << std::endl;
      },
      comm);

  for (const auto &i : sources)
    deallog << i << ' ';
  deallog << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const MPI_Comm comm = MPI_COMM_WORLD;

  test(comm);
}
