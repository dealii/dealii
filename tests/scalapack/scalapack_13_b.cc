// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
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

// test addition of distributed ScaLAPACKMatrices: A = alpha A + beta B^T

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/scalapack.h>

#include <fstream>
#include <iostream>
#include <typeinfo>


template <typename NumberType>
void
test(const unsigned int block_size_i, const unsigned int block_size_j)
{
  MPI_Comm           mpi_communicator(MPI_COMM_WORLD);
  const unsigned int n_mpi_processes(
    Utilities::MPI::n_mpi_processes(mpi_communicator));
  const unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  std::cout << std::setprecision(10);
  ConditionalOStream pcout(std::cout, (this_mpi_process == 0));

  const auto proc_rows =
    static_cast<unsigned int>(std::floor(std::sqrt(n_mpi_processes)));
  const auto proc_columns =
    static_cast<unsigned int>(std::floor(n_mpi_processes / proc_rows));
  // create 2d process grid
  std::shared_ptr<Utilities::MPI::ProcessGrid> grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(mpi_communicator,
                                                  proc_rows,
                                                  proc_columns);
  pcout << "2D process grid: " << grid->get_process_grid_rows() << 'x'
        << grid->get_process_grid_columns() << std::endl
        << std::endl;

  const std::vector<unsigned int> sizes = {{400, 500}};

  FullMatrix<NumberType> full_A(sizes[0], sizes[1]);
  FullMatrix<NumberType> full_B(sizes[1], sizes[0]);
  create_random(full_A);
  create_random(full_B);

  // conditions for block sizes: mb_A=mb_C, nb_B=nb_C, nb_A=mb_B
  const unsigned int mb_A = block_size_i, nb_A = block_size_j;
  const unsigned int mb_B = nb_A, nb_B = mb_A;

  ScaLAPACKMatrix<NumberType> scalapack_A(
    full_A.m(), full_A.n(), grid, mb_A, nb_A);
  ScaLAPACKMatrix<NumberType> scalapack_B(
    full_B.m(), full_B.n(), grid, mb_B, nb_B);
  scalapack_A = full_A;
  scalapack_B = full_B;

  const NumberType alpha = 1.2, beta = -0.7;

  full_A *= alpha;
  FullMatrix<NumberType> full_B_t(sizes[0], sizes[1]);
  full_B_t.copy_transposed(full_B);
  full_A.add(beta, full_B_t);

  scalapack_A.add(scalapack_B, alpha, beta, true);
  FullMatrix<NumberType> tmp_full_A(scalapack_A.m(), scalapack_A.n());
  scalapack_A.copy_to(tmp_full_A);

  pcout << "   computing A = alpha A + beta B^T with"
        << " A in R^(" << scalapack_A.m() << 'x' << scalapack_A.n() << ") and"
        << " B in R^(" << scalapack_B.m() << 'x' << scalapack_B.n() << ')'
        << std::endl;
  pcout << "   norms: " << tmp_full_A.frobenius_norm() << " & "
        << full_A.frobenius_norm() << "  for " << typeid(NumberType).name()
        << std::endl
        << std::endl;
  pcout << std::endl;
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
      test<double>(s, b);
}
