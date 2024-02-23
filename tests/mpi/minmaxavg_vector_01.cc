// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check Utilities::MPI::min_max_avg() for std::vector and ArrayView

#include <deal.II/base/utilities.h>

#include "../tests.h"

void
print_it(const Utilities::MPI::MinMaxAvg &result)
{
  deallog << "sum: " << result.sum << " avg: " << result.avg
          << " min: " << result.min << " @" << result.min_index
          << " max: " << result.max << " @" << result.max_index << std::endl;
}

void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  std::vector<double> values;

  if (myid == 0)
    values = std::vector<double>{1.0, 1.0, 3.0};
  else if (myid == 1)
    values = std::vector<double>{2.0, 3.0, 2.0};
  else if (myid == 2)
    values = std::vector<double>{3.0, 2.0, 1.0};
  else if (myid == 3)
    values = std::vector<double>{3.0, 3.0, 0.0};

  std::vector<Utilities::MPI::MinMaxAvg> results(values.size());

  // test ArrayView
  Utilities::MPI::min_max_avg(values, results, MPI_COMM_WORLD);
  for (const auto &result : results)
    print_it(result);

  // test std::vector
  for (const auto &result : Utilities::MPI::min_max_avg(values, MPI_COMM_WORLD))
    print_it(result);
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test();
}
