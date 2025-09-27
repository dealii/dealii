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



// Test Utilities::MPI::broadcast() for built-in types.

#include <deal.II/base/mpi.h>

#include <iomanip>

#include "../tests.h"

void
check(const int int_value)
{
  const auto my_proc = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  int buffer;
  if (my_proc == 0)
    {
      buffer = int_value;
    }

  const auto result = Utilities::MPI::broadcast(MPI_COMM_WORLD, buffer);

  deallog << result << " ";
  deallog << std::endl;
}


void
check(const double double_value)
{
  const auto my_proc = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  double buffer;
  if (my_proc == 0)
    {
      buffer = double_value;
    }

  const auto result = Utilities::MPI::broadcast(MPI_COMM_WORLD, buffer);

  deallog << std::setprecision(4) << result << " ";
  deallog << std::endl;
}


void
check(const std::complex<double> &complex_value)
{
  const auto my_proc = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  std::complex<double> buffer;
  if (my_proc == 0)
    {
      buffer = complex_value;
    }

  const auto result = Utilities::MPI::broadcast(MPI_COMM_WORLD, buffer);

  deallog << std::setprecision(2) << result << " ";
  deallog << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  const unsigned int root = 0;
  const bool         on_root =
    (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == root);

  // Test broadcast with int, double, and std::complex data
  // types. Only provide the right value on the root, and some other
  // value on non-root processes. We should only see the root value in
  // the output.
  check(on_root ? 23 : 12);
  check(on_root ? 14.32 : std::numeric_limits<double>::max());
  check(on_root ? std::complex<double>(1., 2.) :
                  std::complex<double>(1., 222222.));
}
