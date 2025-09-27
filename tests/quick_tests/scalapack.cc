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

// adaptation of /scalapack/scalapack_02.cc to a quick test

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/scalapack.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iostream>



template <typename FullMatrix>
void
create_spd(FullMatrix &A)
{
  const unsigned int size = A.n();
  Assert(size == A.m(), ExcDimensionMismatch(size, A.m()));

  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = i; j < size; ++j)
      {
        const double val = random_value<typename FullMatrix::value_type>();
        Assert(val >= 0. && val <= 1., ExcInternalError());
        if (i == j)
          // since A(i,j) < 1 and
          // a symmetric diagonally dominant matrix is SPD
          A(i, j) = val + size;
        else
          {
            A(i, j) = val;
            A(j, i) = val;
          }
      }
}



template <typename NumberType>
void
test(const unsigned int size, const unsigned int block_size)
{
  MPI_Comm           mpi_communicator(MPI_COMM_WORLD);
  const unsigned int n_mpi_processes(
    Utilities::MPI::n_mpi_processes(mpi_communicator));
  const unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  ConditionalOStream pcout(std::cout, (this_mpi_process == 0));

  // Create SPD matrices of requested size:
  FullMatrix<NumberType> full_in(size), inverse(size), full_out(size),
    diff(size);

  std::shared_ptr<Utilities::MPI::ProcessGrid> grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(
      mpi_communicator, size, size, block_size, block_size);
  ScaLAPACKMatrix<NumberType> scalapack_matrix(
    size, grid, block_size, LAPACKSupport::Property::symmetric);

  pcout << size << ' ' << block_size << ' ' << grid->get_process_grid_rows()
        << ' ' << grid->get_process_grid_columns() << std::endl;

  create_spd(full_in);

  // invert via Lapack
  inverse.cholesky(full_in);

  // invert via ScaLAPACK
  scalapack_matrix = full_in;
  scalapack_matrix.compute_cholesky_factorization();
  scalapack_matrix.copy_to(full_out);

  diff = 0;
  diff.add(1., inverse);
  diff.add(-1., full_out);

  const double error = diff.frobenius_norm();

  if (error > 1e-10 && this_mpi_process == 0)
    {
      std::cout << "Error!" << std::endl << "expected to have:" << std::endl;
      inverse.print_formatted(std::cout);
      std::cout << "but got:" << std::endl;
      full_out.print_formatted(std::cout);
      std::cout << "difference:" << std::endl;
      diff.print_formatted(std::cout);
      AssertThrow(false, dealii::ExcInternalError());
    }
  else
    pcout << "Ok" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  test<double>(320, 64);
}
