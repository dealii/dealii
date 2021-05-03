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

// test Cholesky decomposition in ScaLAPACK vs FullMatrix

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
  FullMatrix<NumberType> full_in(size), inverse(size), full_out(size),
    diff(size);

  std::shared_ptr<Utilities::MPI::ProcessGrid> grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(
      mpi_communicator, size, size, block_size, block_size);
  ScaLAPACKMatrix<NumberType> scalapack_matrix(
    size, grid, block_size, LAPACKSupport::Property::symmetric);

  pcout << size << " " << block_size << " " << grid->get_process_grid_rows()
        << " " << grid->get_process_grid_columns() << std::endl;

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

  const NumberType error = diff.frobenius_norm() / inverse.frobenius_norm();

  if (error > tol && this_mpi_process == 0)
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

  const std::vector<unsigned int> sizes  = {{32, 64, 120, 320, 640}};
  const std::vector<unsigned int> blocks = {{32, 64}};

  const double tol_double = 1e-10;
  const float  tol_float  = 1e-5;

  for (const auto &s : sizes)
    for (const auto &b : blocks)
      if (b <= s)
        test<float>(s, b, tol_float);

  for (const auto &s : sizes)
    for (const auto &b : blocks)
      if (b <= s)
        test<double>(s, b, tol_double);
}
