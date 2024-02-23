// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test Utilities::MPI::all_reduce().

#include <deal.II/base/mpi.h>

#include "../tests.h"


template <typename T>
void
check(const std::function<std::vector<T>(const std::vector<T> &,
                                         const std::vector<T> &)> &fu)
{
  const auto result = Utilities::MPI::all_reduce(
    std::vector<T>{Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)},
    MPI_COMM_WORLD,
    fu);

  for (const auto r : result)
    deallog << r << ' ';
  deallog << std::endl;

  for (unsigned int rank = 0;
       rank < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
       ++rank)
    {
      const auto result = Utilities::MPI::reduce(
        std::vector<T>{Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)},
        MPI_COMM_WORLD,
        fu,
        rank);

      for (const auto r : result)
        deallog << r << ' ';
      deallog << std::endl;
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  // compute min
  check<unsigned int>([](const auto &a, const auto &b) {
    return std::vector<unsigned int>{std::min(a[0], b[0])};
  });

  // compute max
  check<unsigned int>([](const auto &a, const auto &b) {
    return std::vector<unsigned int>{std::max(a[0], b[0])};
  });

  // compute sum
  check<unsigned int>([](const auto &a, const auto &b) {
    return std::vector<unsigned int>{a[0] + b[0]};
  });

  // perform gather all
  check<unsigned int>([](const auto &a, const auto &b) {
    auto result = a;
    result.insert(result.end(), b.begin(), b.end());
    return result;
  });
}
