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

// test copying of distributed ScaLAPACKMatrices using ScaLAPACK routine
// p_gemr2d

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/scalapack.h>

#include <fstream>
#include <iostream>


template <typename NumberType>
void
test(const unsigned int block_size_i, const unsigned int block_size_j)
{
  MPI_Comm           mpi_communicator(MPI_COMM_WORLD);
  const unsigned int n_mpi_processes(
    Utilities::MPI::n_mpi_processes(mpi_communicator));
  const unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  ConditionalOStream pcout(std::cout, (this_mpi_process == 0));

  const unsigned int size = 500;
  // create FullMatrix and fill it
  FullMatrix<NumberType> full(size);
  unsigned int           count = 0;
  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j, ++count)
      full(i, j) = count;

  // create 2d process grid
  std::shared_ptr<Utilities::MPI::ProcessGrid> grid_2d =
    std::make_shared<Utilities::MPI::ProcessGrid>(
      mpi_communicator, size, size, block_size_i, block_size_i);
  // create 1d process grid
  std::shared_ptr<Utilities::MPI::ProcessGrid> grid_1d =
    std::make_shared<Utilities::MPI::ProcessGrid>(mpi_communicator,
                                                  n_mpi_processes,
                                                  1);
  // create process grid containing one process
  std::shared_ptr<Utilities::MPI::ProcessGrid> grid_single =
    std::make_shared<Utilities::MPI::ProcessGrid>(mpi_communicator, 1, 1);

  ScaLAPACKMatrix<NumberType> scalapack_matrix_2d(
    size, size, grid_2d, block_size_i, block_size_i);
  ScaLAPACKMatrix<NumberType> scalapack_matrix_1d(
    size, size, grid_1d, block_size_j, block_size_j);
  ScaLAPACKMatrix<NumberType> scalapack_matrix_single(
    size, size, grid_single, block_size_i, block_size_j);
  ScaLAPACKMatrix<NumberType> scalapack_matrix_source(
    size, size, grid_single, block_size_j, block_size_i);

  pcout << "2D grid matrix: dim=" << scalapack_matrix_2d.m() << 'x'
        << scalapack_matrix_2d.n() << ";  blocks=" << block_size_i << 'x'
        << block_size_i << ";  grid=" << grid_2d->get_process_grid_rows() << 'x'
        << grid_2d->get_process_grid_columns() << std::endl;

  pcout << "1D grid matrix: " << scalapack_matrix_1d.m() << 'x'
        << scalapack_matrix_1d.n() << ";  blocks=" << block_size_j << 'x'
        << block_size_j << ";  grid=" << grid_1d->get_process_grid_rows() << 'x'
        << grid_1d->get_process_grid_columns() << std::endl;

  pcout << "single process matrix: " << scalapack_matrix_single.m() << 'x'
        << scalapack_matrix_single.n() << ";  blocks=" << block_size_i << 'x'
        << block_size_j << ";  grid=" << grid_single->get_process_grid_rows()
        << 'x' << grid_single->get_process_grid_columns() << std::endl
        << std::endl;

  scalapack_matrix_source = full;

  scalapack_matrix_source.copy_to(scalapack_matrix_2d);
  scalapack_matrix_source.copy_to(scalapack_matrix_1d);
  scalapack_matrix_source.copy_to(scalapack_matrix_single);

  FullMatrix<NumberType> test_2d(size), test_1d(size), test_one(size);
  scalapack_matrix_2d.copy_to(test_2d);
  scalapack_matrix_1d.copy_to(test_1d);
  scalapack_matrix_single.copy_to(test_one);
  test_2d.add(-1, full);
  test_1d.add(-1, full);
  test_one.add(-1, full);

  AssertThrow(test_2d.frobenius_norm() < 1e-12, ExcInternalError());
  AssertThrow(test_1d.frobenius_norm() < 1e-12, ExcInternalError());
  AssertThrow(test_one.frobenius_norm() < 1e-12, ExcInternalError());
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  const std::vector<unsigned int> blocks_i = {{16, 32, 64}};
  const std::vector<unsigned int> blocks_j = {{16, 32, 64}};

  for (const auto &s : blocks_i)
    for (const auto &b : blocks_j)
      test<float>(s, b);

  for (const auto &s : blocks_i)
    for (const auto &b : blocks_j)
      test<double>(s, b);
}
