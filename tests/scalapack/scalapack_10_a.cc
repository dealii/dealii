// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include "../tests.h"

#include "../lapack/create_matrix.h"

// A test for saving and loading of distributed ScaLAPACKMatrices.
// By using a 1x1 grid a call to save_serial() and load_serial() is mimicked.
// The ScaLAPACKMatrix is saved in parallel/serial and afterwards loaded
// in serial/parallel and the content is compared.
// A former bug (removed in #7044) in save_parallel and load_parallel caused
// an inconsistency so that matrices saved in serial could not be read
// correctly in parallel. This test checks for that interoperability.

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
  const std::string filename_serial("scalapack_10_a_s_test.h5");
  const std::string filename_parallel("scalapack_10_a_p_test.h5");

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

  // create process grid containing only one process
  std::shared_ptr<Utilities::MPI::ProcessGrid> grid_one =
    std::make_shared<Utilities::MPI::ProcessGrid>(mpi_communicator, 1, 1);

  ScaLAPACKMatrix<NumberType> matrix(size, size, grid, block_size, block_size);
  ScaLAPACKMatrix<NumberType> matrix_one(
    size, size, grid, block_size, block_size);

  matrix     = full;
  matrix_one = full;

  matrix.save(filename_parallel);
  matrix_one.save(filename_serial);

  matrix.load(filename_serial);
  matrix_one.load(filename_parallel);

  std::remove(filename_serial.c_str());
  std::remove(filename_parallel.c_str());

  FullMatrix<NumberType> copy_I(size), copy_II(size);
  matrix.copy_to(copy_I);
  matrix_one.copy_to(copy_II);
  copy_I.add(-1, full);
  copy_II.add(-1, full);

  pcout << size << ' ' << block_size << std::endl;

  if (copy_I.frobenius_norm() > 1e-12)
    pcout << "norm of difference I: " << copy_I.frobenius_norm() << std::endl;
  AssertThrow(copy_I.frobenius_norm() < 1e-12, ExcInternalError());

  if (copy_II.frobenius_norm() > 1e-12)
    pcout << "norm of difference II: " << copy_II.frobenius_norm() << std::endl;
  AssertThrow(copy_II.frobenius_norm() < 1e-12, ExcInternalError());
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

  const std::vector<unsigned int> sizes       = {{200, 300}};
  const std::vector<unsigned int> block_sizes = {{16, 32}};

  for (const auto &s : sizes)
    for (const auto &b : block_sizes)
      test<double>(s, b);

  for (const auto &s : sizes)
    for (const auto &b : block_sizes)
      test<float>(s, b);
}
