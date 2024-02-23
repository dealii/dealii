// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
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

// test least_squares(ScaLAPACKMatrix<NumberType>&,const bool)

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/process_grid.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/scalapack.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>


template <typename NumberType>
void
test(const unsigned int block_size, const NumberType tol)
{
  MPI_Comm           mpi_communicator(MPI_COMM_WORLD);
  const unsigned int n_mpi_processes(
    Utilities::MPI::n_mpi_processes(mpi_communicator));
  const unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  ConditionalOStream pcout(std::cout, (this_mpi_process == 0));

  std::shared_ptr<Utilities::MPI::ProcessGrid> grid_2d =
    std::make_shared<Utilities::MPI::ProcessGrid>(
      mpi_communicator, 4, 3, block_size, block_size);

  // examples from
  // https://www.ibm.com/support/knowledgecenter/en/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lgels.htm
  FullMatrix<NumberType> full_A_I(4, 3), full_B_I(4, 5), full_X_I(4, 5);

  // FIXME: Add more tests for different cases!!!
  pcout << "Solving least squares problem ||B - A*X||" << std::endl;

  full_A_I(0, 0) = 1.;
  full_A_I(0, 1) = -2.;
  full_A_I(0, 2) = -1.;
  full_A_I(1, 0) = 2.;
  full_A_I(1, 1) = 0.;
  full_A_I(1, 2) = 1.;
  full_A_I(2, 0) = 2.;
  full_A_I(2, 1) = -4.;
  full_A_I(2, 2) = 2.;
  full_A_I(3, 0) = 4.;
  full_A_I(3, 1) = 0.;
  full_A_I(3, 2) = 0.;

  full_B_I(0, 0) = -1.;
  full_B_I(0, 1) = -2.;
  full_B_I(0, 2) = -7.;
  full_B_I(0, 3) = 0.;
  full_B_I(0, 4) = -5.;
  full_B_I(1, 0) = 1.;
  full_B_I(1, 1) = 3.;
  full_B_I(1, 2) = 4.;
  full_B_I(1, 3) = 3.;
  full_B_I(1, 4) = 5.;
  full_B_I(2, 0) = 1.;
  full_B_I(2, 1) = 0.;
  full_B_I(2, 2) = 4.;
  full_B_I(2, 3) = 2.;
  full_B_I(2, 4) = 2.;
  full_B_I(3, 0) = -2.;
  full_B_I(3, 1) = 4.;
  full_B_I(3, 2) = 4.;
  full_B_I(3, 3) = 0.;
  full_B_I(3, 4) = 4.;

  full_X_I(0, 0) = -0.4;
  full_X_I(0, 1) = 1.;
  full_X_I(0, 2) = 0.8;
  full_X_I(0, 3) = 0.2;
  full_X_I(0, 4) = 1.;
  full_X_I(1, 0) = 0.;
  full_X_I(1, 1) = 1.;
  full_X_I(1, 2) = 1.5;
  full_X_I(1, 3) = 0.;
  full_X_I(1, 4) = 1.5;
  full_X_I(2, 0) = 1.;
  full_X_I(2, 1) = 1.;
  full_X_I(2, 2) = 4.;
  full_X_I(2, 3) = 1.;
  full_X_I(2, 4) = 3.;
  full_X_I(3, 0) = -1.;
  full_X_I(3, 1) = 0.;
  full_X_I(3, 2) = 2.;
  full_X_I(3, 3) = -2.;
  full_X_I(3, 4) = 0.;

  // compute eigenpairs of s.p.d matrix
  ScaLAPACKMatrix<NumberType> scalapack_A(
    4, 3, grid_2d, block_size, block_size);
  ScaLAPACKMatrix<NumberType> scalapack_B(
    4, 5, grid_2d, block_size, block_size);
  scalapack_A.set_property(LAPACKSupport::Property::general);
  scalapack_A = full_A_I;
  scalapack_B = full_B_I;
  scalapack_A.least_squares(scalapack_B, false);
  FullMatrix<NumberType> result(4, 5);
  scalapack_B.copy_to(result);

  result.add(-1, full_X_I);
  AssertThrow(result.frobenius_norm() < tol,
              ExcMessage("solution deviates from reference"));
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  const std::vector<unsigned int> blocks     = {{1, 2}};
  const double                    tol_double = 1e-10;
  const float                     tol_float  = 1e-5;

  for (const auto &b : blocks)
    test<double>(b, tol_double);

  for (const auto &b : blocks)
    test<float>(b, tol_float);
}
