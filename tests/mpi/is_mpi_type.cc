// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test the Utilities::MPI::is_mpi_type template variable.

#include <deal.II/base/mpi.h>

#include "../tests.h"


void
test()
{
  deallog << std::boolalpha;

  // Verify that the following types are all supported:
  deallog << Utilities::MPI::is_mpi_type<char> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<signed short> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<signed int> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<signed long> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<signed long long> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<signed char> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<unsigned char> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<unsigned short> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<unsigned int> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<unsigned long int> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<unsigned long long> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<float> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<double> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<long double> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<bool> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<std::complex<float>> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<std::complex<double>> << std::endl;
  deallog
    << Utilities::MPI::is_mpi_type<std::complex<long double>> << std::endl;
  deallog << Utilities::MPI::is_mpi_type<wchar_t> << std::endl;


  // Then also check a non-native type:
  struct X
  {};
  deallog << Utilities::MPI::is_mpi_type<X> << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test();
}
