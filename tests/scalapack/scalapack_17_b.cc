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

// test ScaLAPACKMatrix::copy_to (LAPACKFullMatrix<NumberType> &matrix)

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

  unsigned int count = 0;
  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j, ++count)
      full(i, j) = count;

  // create 2d process grid
  std::shared_ptr<Utilities::MPI::ProcessGrid> grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(
      mpi_communicator, size, size, block_size_i, block_size_i);
  ScaLAPACKMatrix<NumberType> scalapack_matrix(
    size, size, grid, block_size_i, block_size_i);

  pcout << "2D grid matrix: dim=" << scalapack_matrix.m() << 'x'
        << scalapack_matrix.n() << ";  blocks=" << block_size_i << 'x'
        << block_size_i << ";  grid=" << grid->get_process_grid_rows() << 'x'
        << grid->get_process_grid_columns() << std::endl
        << std::endl;

  scalapack_matrix = full;
  LAPACKFullMatrix<NumberType> comparison;

  const bool       process_is_in_grid = grid->is_process_active();
  const int        entry = process_is_in_grid ? this_mpi_process : -1;
  std::vector<int> valid_process_ranks =
    Utilities::MPI::gather(mpi_communicator, entry, 0);
  int random_rank = 0;

  if (this_mpi_process == 0)
    {
      std::vector<int> valid_ranks;

      for (unsigned int i = 0; i < valid_process_ranks.size(); ++i)
        if (valid_process_ranks[i] >= 0)
          valid_ranks.push_back(valid_process_ranks[i]);

      const int entry = random_value<unsigned int>(0, valid_ranks.size() - 1);
      random_rank     = valid_ranks[entry];
    }
  MPI_Bcast(&random_rank, 1, MPI_INT, 0, mpi_communicator);
  const unsigned int rank = random_rank;

  if (this_mpi_process == rank)
    comparison.reinit(size);

  scalapack_matrix.copy_to(comparison, rank);

  if (this_mpi_process == rank)
    {
      NumberType diff = 0;
      for (unsigned int i = 0; i < size; ++i)
        for (unsigned int j = 0; j < size; ++j)
          diff +=
            (full(i, j) - comparison(i, j)) * (full(i, j) - comparison(i, j));
      diff = sqrt(diff);

      AssertThrow(diff < 1e-12, ExcInternalError());
    }
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
