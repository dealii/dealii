// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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
