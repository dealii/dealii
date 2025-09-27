// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check Utilities::MPI::sum() for FullMatrix objects

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

template <typename NumberType>
void
test(const unsigned int m = 13, const unsigned int n = 5)
{
  Assert(Utilities::MPI::job_supports_mpi(), ExcInternalError());

  FullMatrix<NumberType> full_matrix(m, n);
  {
    unsigned int index = 0;
    for (unsigned int i = 0; i < full_matrix.m(); ++i)
      for (unsigned int j = 0; j < full_matrix.n(); ++j)
        full_matrix(i, j) = index++;
  }

  FullMatrix<NumberType> full_matrix_original(m, n);
  full_matrix_original = full_matrix;

  // inplace
  Utilities::MPI::sum(full_matrix, MPI_COMM_WORLD, full_matrix);

  const unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  for (unsigned int i = 0; i < full_matrix.m(); ++i)
    for (unsigned int j = 0; j < full_matrix.n(); ++j)
      Assert(full_matrix(i, j) == full_matrix_original(i, j) * double(numprocs),
             ExcInternalError());

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "Ok" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();
      deallog.push("float");
      test<float>();
      deallog.pop();
      deallog.push("double");
      test<double>();
      deallog.pop();
#ifdef DEAL_II_WITH_COMPLEX_VALUES
      deallog.push("complex<double>");
      test<std::complex<double>>();
#endif
      deallog.pop();
    }
  else
    {
      test<float>();
      test<double>();
#ifdef DEAL_II_WITH_COMPLEX_VALUES
      test<std::complex<double>>();
#endif
    }
}
