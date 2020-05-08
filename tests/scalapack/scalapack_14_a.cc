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

// test scaling of rows of distributed ScaLAPACKMatrices

#include <deal.II/base/array_view.h>
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
  const std::vector<unsigned int>              sizes = {{400, 500}};
  std::shared_ptr<Utilities::MPI::ProcessGrid> grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(
      mpi_communicator, sizes[0], sizes[1], block_size_i, block_size_i);
  pcout << "2D process grid: " << grid->get_process_grid_rows() << "x"
        << grid->get_process_grid_columns() << std::endl
        << std::endl;

  // test scaling of rows
  FullMatrix<NumberType> full_A(sizes[0], sizes[1]);
  create_random(full_A);

  std::vector<NumberType> scaling_factors(full_A.m());
  for (unsigned int i = 0; i < scaling_factors.size(); ++i)
    scaling_factors[i] = std::sqrt(i + 1);

  ScaLAPACKMatrix<NumberType> scalapack_A(
    full_A.m(), full_A.n(), grid, block_size_i, block_size_j);
  scalapack_A = full_A;
  const ArrayView<NumberType> view_rows(scaling_factors);
  scalapack_A.scale_rows(view_rows);
  FullMatrix<NumberType> tmp_full_A(scalapack_A.m(), scalapack_A.n());
  scalapack_A.copy_to(tmp_full_A);

  for (unsigned int i = 0; i < full_A.m(); ++i)
    for (unsigned int j = 0; j < full_A.n(); ++j)
      full_A(i, j) *= scaling_factors[i];

  pcout << "   Row scaling for"
        << " A in R^(" << scalapack_A.m() << "x" << scalapack_A.n() << ")"
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
