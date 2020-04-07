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


// Test the function compute_set_union() for set::vector and std::set.

#include <deal.II/base/mpi.h>

#include "../tests.h"

using namespace dealii;

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
      deallog << i << " ";
    deallog << std::endl;
  }

  // test function for set
  {
    const auto result = Utilities::MPI::compute_set_union(
      std::set<unsigned int>(vector.begin(), vector.end()), comm);

    for (auto i : result)
      deallog << i << " ";
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
