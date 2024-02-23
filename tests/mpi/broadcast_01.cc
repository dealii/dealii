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



// Test Utilities::MPI::broadcast().

#include <deal.II/base/mpi.h>

#include "../tests.h"

struct TestObject
{
  std::vector<int> int_vector;
  std::string      string;

  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int /*version*/)
  {
    ar &int_vector;
    ar &string;
  }

  friend std::ostream &
  operator<<(std::ostream &os, const TestObject &test_object);
};

std::ostream &
operator<<(std::ostream &os, const TestObject &test_object)
{
  for (const auto &i : test_object.int_vector)
    os << i << ' ';

  os << test_object.string;

  return os;
}

void
check(const std::vector<int> &int_vector, const std::string &string)
{
  const auto my_proc = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  TestObject test_object;

  if (my_proc == 0)
    {
      test_object.int_vector = int_vector;
      test_object.string     = string;
    }

  const auto result = Utilities::MPI::broadcast(MPI_COMM_WORLD, test_object);

  deallog << result << ' ';
  deallog << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  // broadcast test
  check({1, 2, 3}, "test");

  // broadcast test
  check({-1, 0, 1}, "i love democracy");
}
