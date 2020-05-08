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

// test norms of ScaLAPACK vs FullMatrix for general matrices

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
test(const unsigned int size,
     const unsigned int block_size,
     const NumberType   tol)
{
  MPI_Comm           mpi_communicator(MPI_COMM_WORLD);
  const unsigned int n_mpi_processes(
    Utilities::MPI::n_mpi_processes(mpi_communicator));
  const unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  ConditionalOStream pcout(std::cout, (this_mpi_process == 0));

  // Create SPD matrices of requested size:
  FullMatrix<NumberType> full(size);

  std::shared_ptr<Utilities::MPI::ProcessGrid> grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(
      mpi_communicator, size, size, block_size, block_size);
  ScaLAPACKMatrix<NumberType> scalapack(
    size, size, grid, block_size, block_size, LAPACKSupport::Property::general);

  pcout << size << " " << block_size << " " << grid->get_process_grid_rows()
        << " " << grid->get_process_grid_columns() << std::endl;

  create_random(full);
  scalapack = full;

  const NumberType l1        = full.l1_norm();
  const NumberType linfty    = full.linfty_norm();
  const NumberType frobenius = full.frobenius_norm();

  // local result on this core:
  const NumberType s_l1        = scalapack.l1_norm();
  const NumberType s_linfty    = scalapack.linfty_norm();
  const NumberType s_frobenius = scalapack.frobenius_norm();

  // make sure we have the same result on all cores, do average:
  const NumberType as_l1 =
    dealii::Utilities::MPI::sum(s_l1, mpi_communicator) / n_mpi_processes;
  const NumberType as_linfty =
    dealii::Utilities::MPI::sum(s_linfty, mpi_communicator) / n_mpi_processes;
  const NumberType as_frobenius =
    dealii::Utilities::MPI::sum(s_frobenius, mpi_communicator) /
    n_mpi_processes;

  pcout << l1 << " " << s_l1 << " " << as_l1 << std::endl
        << linfty << " " << s_linfty << " " << as_linfty << std::endl
        << frobenius << " " << s_frobenius << " " << as_frobenius << std::endl;

  AssertThrow(std::abs(l1 - as_l1) < tol * std::abs(l1),
              dealii::ExcInternalError());
  AssertThrow(std::abs(linfty - as_linfty) < tol * std::abs(linfty),
              dealii::ExcInternalError());
  AssertThrow(std::abs(frobenius - as_frobenius) < tol * std::abs(frobenius),
              dealii::ExcInternalError());
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  const std::vector<unsigned int> sizes  = {{32, 64, 120, 320, 640}};
  const std::vector<unsigned int> blocks = {{32, 64}};

  const double tol_double = 1e-10;
  const float  tol_float  = 1e-5;

  /*for (const auto &s : sizes)
    for (const auto &b : blocks)
      if (b <= s)
        test<float>(s,b,tol_float);*/

  for (const auto &s : sizes)
    for (const auto &b : blocks)
      if (b <= s)
        test<double>(s, b, tol_double);
}
