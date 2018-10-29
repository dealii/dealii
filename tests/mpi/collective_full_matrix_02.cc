// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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



// check Utilities::MPI::sum() for LAPACKFullMatrix objects

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

template <typename NumberType>
void
test(const unsigned int m = 13, const unsigned int n = 5)
{
  Assert(Utilities::MPI::job_supports_mpi(), ExcInternalError());

  LAPACKFullMatrix<NumberType> full_matrix(m, n);
  {
    unsigned int index = 0;
    for (unsigned int i = 0; i < full_matrix.m(); ++i)
      for (unsigned int j = 0; j < full_matrix.n(); ++j)
        full_matrix(i, j) = index++;
  }

  LAPACKFullMatrix<NumberType> full_matrix_original(m, n);
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
    }
  else
    {
      test<float>();
      test<double>();
    }
}
