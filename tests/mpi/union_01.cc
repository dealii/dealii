// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test the function compute_set_union() for set::vector and std::set.

#include <deal.II/base/mpi.h>

#include "../tests.h"


void
test(const MPI_Comm comm)
{
  std::vector<unsigned int> vector;

  if (Utilities::MPI::this_mpi_process(comm) == 0)
    vector = std::vector<unsigned int>{2, 5, 7};
  else
    vector = std::vector<unsigned int>{5, 8, 9, 10};

  // test function for vector
  {
    const auto result = Utilities::MPI::compute_set_union(vector, comm);

    for (auto i : result)
      deallog << i << ' ';
    deallog << std::endl;
  }

  // test function for set
  {
    const auto result = Utilities::MPI::compute_set_union(
      std::set<unsigned int>(vector.begin(), vector.end()), comm);

    for (auto i : result)
      deallog << i << ' ';
    deallog << std::endl;
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
