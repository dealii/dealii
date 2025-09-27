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

// test multiplication of distributed ScaLAPACKMatrices: C = alpha A^T*B + beta
// C

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
test()
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

  const std::vector<unsigned int> sizes = {{300, 400, 500}};

  FullMatrix<NumberType> full_A(sizes[2], sizes[0]);
  FullMatrix<NumberType> full_B(sizes[2], sizes[1]);
  FullMatrix<NumberType> full_C(sizes[0], sizes[1]);
  create_random(full_A);
  create_random(full_B);
  create_random(full_C);

  // conditions for block sizes: nb_A=mb_C, nb_B=nb_C, mb_A=mb_B
  const unsigned int mb_A = 32, nb_A = 64, nb_B = 16;
  const unsigned int mb_B = mb_A, mb_C = nb_A;
  const unsigned int nb_C = nb_B;

  ScaLAPACKMatrix<NumberType> scalapack_A(
    full_A.m(), full_A.n(), grid, mb_A, nb_A);
  ScaLAPACKMatrix<NumberType> scalapack_B(
    full_B.m(), full_B.n(), grid, mb_B, nb_B);
  ScaLAPACKMatrix<NumberType> scalapack_C(
    full_C.m(), full_C.n(), grid, mb_C, nb_C);
  scalapack_A = full_A;
  scalapack_B = full_B;
  scalapack_C = full_C;

  const NumberType b = 1.4, c = 0.1;

  full_A *= b;
  full_C *= c;
  full_A.Tmmult(full_C, full_B, true);

  scalapack_A.mult(b, scalapack_B, c, scalapack_C, true, false);
  FullMatrix<NumberType> tmp_full_C(full_C.m(), full_C.n());
  scalapack_C.copy_to(tmp_full_C);

  pcout << "   computing C = b A^T * B + c C with"
        << " A in R^(" << scalapack_A.m() << 'x' << scalapack_A.n() << "),"
        << " B in R^(" << scalapack_B.m() << 'x' << scalapack_B.n() << ") and"
        << " C in R^(" << scalapack_C.m() << 'x' << scalapack_C.n() << ')'
        << std::endl;
  pcout << "   norms: " << tmp_full_C.frobenius_norm() << " & "
        << full_C.frobenius_norm() << "  for " << typeid(NumberType).name()
        << std::endl
        << std::endl;
  pcout << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  test<double>();
}
