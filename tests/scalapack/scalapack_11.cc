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

// test copying submatrices of distributed ScaLAPACKMatrices using ScaLAPACK
// routine p_gemr2d

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
  std::shared_ptr<Utilities::MPI::ProcessGrid> grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(
      mpi_communicator, size, size, block_size_i, block_size_i);

  ScaLAPACKMatrix<NumberType> scalapack_matrix(
    size, size, grid, block_size_i, block_size_i);

  pcout << "2D grid matrix: dim=" << scalapack_matrix.m() << "x"
        << scalapack_matrix.n() << ";  blocks=" << block_size_i << "x"
        << block_size_i << ";  grid=" << grid->get_process_grid_rows() << "x"
        << grid->get_process_grid_columns() << std::endl
        << std::endl;

  scalapack_matrix                               = full;
  unsigned int                          sub_size = 100;
  std::pair<unsigned int, unsigned int> offset_A = std::make_pair(49, 99);
  std::pair<unsigned int, unsigned int> offset_B = std::make_pair(4, 7);
  std::pair<unsigned int, unsigned int> submatrix_size =
    std::make_pair(sub_size, sub_size);
  ScaLAPACKMatrix<NumberType> scalapack_matrix_dest(sub_size + offset_B.first,
                                                    sub_size + offset_B.second,
                                                    grid,
                                                    block_size_j,
                                                    block_size_i);
  scalapack_matrix.copy_to(scalapack_matrix_dest,
                           offset_A,
                           offset_B,
                           submatrix_size);
  FullMatrix<NumberType> dest(sub_size + offset_B.first,
                              sub_size + offset_B.second);
  scalapack_matrix_dest.copy_to(dest);

  for (unsigned int i = 0; i < sub_size; ++i)
    for (unsigned int j = 0; j < sub_size; ++j)
      dest(i + offset_B.first, j + offset_B.second) -=
        full(offset_A.first + i, offset_A.second + j);
  AssertThrow(dest.frobenius_norm() < 1e-12, ExcInternalError());
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
