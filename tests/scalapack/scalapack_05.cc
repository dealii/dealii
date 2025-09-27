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

// test reciprocal_condition_number()

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/scalapack.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iostream>

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
  pcout << std::setprecision(7);

  // Create SPD matrices of requested size:
  FullMatrix<NumberType> full_A(size), inv_A(size);

  std::shared_ptr<Utilities::MPI::ProcessGrid> grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(
      mpi_communicator, size, size, block_size, block_size);
  ScaLAPACKMatrix<NumberType> scalapack_A(size,
                                          grid,
                                          block_size,
                                          LAPACKSupport::Property::symmetric);

  pcout << size << ' ' << block_size << ' ' << grid->get_process_grid_rows()
        << ' ' << grid->get_process_grid_columns() << std::endl;

  create_spd(full_A);
  inv_A.invert(full_A);

  const NumberType l1     = full_A.l1_norm();
  const NumberType inv_l1 = inv_A.l1_norm();

  // Scalapack:
  scalapack_A                   = full_A;
  const NumberType scalapack_l1 = scalapack_A.l1_norm();
  scalapack_A.compute_cholesky_factorization();
  const NumberType rcond =
    scalapack_A.reciprocal_condition_number(scalapack_l1);


  pcout << 1. / (l1 * inv_l1) << ' ' << rcond << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  const std::vector<unsigned int> sizes  = {{32, 64, 120, 320, 640}};
  const std::vector<unsigned int> blocks = {{32, 64}};

  for (const auto &s : sizes)
    for (const auto &b : blocks)
      if (b <= s)
        test<float>(s, b);


  for (const auto &s : sizes)
    for (const auto &b : blocks)
      if (b <= s)
        test<double>(s, b);
}
