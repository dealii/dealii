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

// test saving and loading of distributed ScaLAPACKMatrices with prescribed
// chunk sizes

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
test(const std::pair<unsigned int, unsigned int> &size,
     const unsigned int                           block_size,
     const std::pair<unsigned int, unsigned int> &chunk_size)
{
  const std::string filename("scalapack_10_b_test.h5");

  MPI_Comm           mpi_communicator(MPI_COMM_WORLD);
  const unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));
  ConditionalOStream pcout(std::cout, (this_mpi_process == 0));

  FullMatrix<NumberType> full(size.first, size.second);
  create_random(full);

  // create 2d process grid
  std::shared_ptr<Utilities::MPI::ProcessGrid> grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(
      mpi_communicator, size.first, size.second, block_size, block_size);

  ScaLAPACKMatrix<NumberType> scalapack_matrix(
    size.first, size.second, grid, block_size, block_size);
  ScaLAPACKMatrix<NumberType> scalapack_matrix_copy(
    size.first, size.second, grid, block_size, block_size);

  scalapack_matrix = full;
  scalapack_matrix.save(filename, chunk_size);
  scalapack_matrix_copy.load(filename);

  FullMatrix<NumberType> copy(size.first, size.second);
  scalapack_matrix_copy.copy_to(copy);
  copy.add(-1, full);

  pcout << size.first << 'x' << size.second << " & " << block_size << " & "
        << chunk_size.first << 'x' << chunk_size.second << " & "
        << grid->get_process_grid_rows() << 'x'
        << grid->get_process_grid_columns() << std::endl;
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

  std::vector<std::pair<unsigned int, unsigned int>> sizes;
  sizes.push_back(std::make_pair(300, 250));

  const std::vector<unsigned int> block_sizes = {{32, 64}};

  std::vector<std::pair<unsigned int, unsigned int>> chunk_sizes;
  chunk_sizes.push_back(std::make_pair(25, 25));
  chunk_sizes.push_back(std::make_pair(75, 50));
  chunk_sizes.push_back(std::make_pair(100, 75));

  for (unsigned int i = 0; i < sizes.size(); ++i)
    for (unsigned int j = 0; j < block_sizes.size(); ++j)
      for (unsigned int k = 0; k < chunk_sizes.size(); ++k)
        test<double>(sizes[i], block_sizes[j], chunk_sizes[k]);

  for (unsigned int i = 0; i < sizes.size(); ++i)
    for (unsigned int j = 0; j < block_sizes.size(); ++j)
      for (unsigned int k = 0; k < chunk_sizes.size(); ++k)
        test<float>(sizes[i], block_sizes[j], chunk_sizes[k]);
}
