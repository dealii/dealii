// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

#include "../tests.h"

#include "../lapack/create_matrix.h"

// test saving and loading of distributed ScaLAPACKMatrices

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/scalapack.h>

#include <cstdio>
#include <fstream>
#include <iostream>


template <typename NumberType>
void
test(const unsigned int size, const unsigned int block_size)
{
  const std::string filename("scalapack_10_test.h5");

  MPI_Comm           mpi_communicator(MPI_COMM_WORLD);
  const unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));
  ConditionalOStream pcout(std::cout, (this_mpi_process == 0));

  // create FullMatrix and fill it
  FullMatrix<NumberType> full(size);
  unsigned int           count = 0;
  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j, ++count)
      full(i, j) = count;

  // create 2d process grid
  std::shared_ptr<Utilities::MPI::ProcessGrid> grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(
      mpi_communicator, size, size, block_size, block_size);

  ScaLAPACKMatrix<NumberType> scalapack_matrix(
    size, size, grid, block_size, block_size);
  ScaLAPACKMatrix<NumberType> scalapack_matrix_copy(
    size, size, grid, block_size, block_size);

  scalapack_matrix = full;
  scalapack_matrix.save(filename);
  scalapack_matrix_copy.load(filename);

  FullMatrix<NumberType> copy(size);
  scalapack_matrix_copy.copy_to(copy);
  copy.add(-1, full);

  pcout << size << " " << block_size << std::endl;

  if (copy.frobenius_norm() > 1e-12)
    pcout << "norm of difference: " << copy.frobenius_norm() << std::endl;

  AssertThrow(copy.frobenius_norm() < 1e-12, ExcInternalError());
  std::remove(filename.c_str());
}



int
main(int argc, char **argv)
{
  // tests.h enables floating point exceptions in debug mode, but this test
  // generates an (irrelevant) exception when run with more than one MPI
  // process so disable them again:
#if defined(DEBUG) && defined(DEAL_II_HAVE_FP_EXCEPTIONS)
  {
    const int current_fe_except = fegetexcept();
    fedisableexcept(current_fe_except);
  }
#endif

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  const std::vector<unsigned int> sizes       = {{100, 200, 300}};
  const std::vector<unsigned int> block_sizes = {{1, 16, 32}};

  for (const auto &s : sizes)
    for (const auto &b : block_sizes)
      test<double>(s, b);

  for (const auto &s : sizes)
    for (const auto &b : block_sizes)
      test<float>(s, b);
}
