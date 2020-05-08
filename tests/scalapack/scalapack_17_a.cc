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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include "../tests.h"

#include "../lapack/create_matrix.h"

// test ScaLAPACKMatrix::operator = (const LAPACKFullMatrix<NumberType> &)

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/scalapack.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>


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

  // create LAPACKFullMatrix and fill it
  LAPACKFullMatrix<NumberType> full;

  if (this_mpi_process == rank)
    {
      full.reinit(size);

      unsigned int count = 0;
      for (unsigned int i = 0; i < size; ++i)
        for (unsigned int j = 0; j < size; ++j, ++count)
          full(i, j) = count;
    }

  scalapack_matrix.copy_from(full, rank);
  FullMatrix<NumberType> comparison(size);
  scalapack_matrix.copy_to(comparison);

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
