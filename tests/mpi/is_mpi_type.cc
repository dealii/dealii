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


// Test the Utilities::MPI::is_mpi_type template variable.

#include <deal.II/base/mpi.h>

#include "../tests.h"

using namespace dealii;

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
