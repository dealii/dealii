// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test ConsensusAlgorithms::selection().

#include <deal.II/base/mpi_consensus_algorithms.h>

#include "../tests.h"


void
test(const MPI_Comm comm)
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
