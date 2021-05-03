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

// test pseudo_inverse(const NumberType ratio) for symmetric matrices

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/process_grid.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/scalapack.h>
#include <deal.II/lac/vector.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>


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

  std::shared_ptr<Utilities::MPI::ProcessGrid> grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(
      mpi_communicator, size, size, block_size, block_size);

  pcout << size << " " << block_size << " " << grid->get_process_grid_rows()
        << " " << grid->get_process_grid_columns() << std::endl;

  // Create SPD matrices of requested size:
  FullMatrix<NumberType> full_A(size), singular_A(size);
  create_spd(full_A);
  singular_A = full_A;
  // Setting last row and column of the matrix to zero to make the matrix
  // singular.
  for (unsigned int i = 0; i < size; ++i)
    {
      singular_A(i, size - 1) = 0;
      singular_A(size - 1, i) = 0;
    }
  const NumberType ratio = 1e-8;

  ScaLAPACKMatrix<NumberType> scalapack_A_1(size, grid, block_size);
  scalapack_A_1 = singular_A;
  scalapack_A_1.set_property(LAPACKSupport::Property::general);
  const unsigned int rank_1 = scalapack_A_1.pseudoinverse(ratio);

  ScaLAPACKMatrix<NumberType> inverse_A(size, grid, block_size);
  inverse_A = full_A;
  inverse_A.set_property(LAPACKSupport::Property::symmetric);
  inverse_A.invert();

  ScaLAPACKMatrix<NumberType> scalapack_A_2(size, grid, block_size);
  scalapack_A_2             = full_A;
  const unsigned int rank_2 = scalapack_A_2.pseudoinverse(ratio);

  FullMatrix<NumberType> full_inverse(size);
  inverse_A.copy_to(full_inverse);
  FullMatrix<NumberType> full_pseudo_inverse(size);
  scalapack_A_2.copy_to(full_pseudo_inverse);

  if (this_mpi_process == 0)
    {
      std::cout << "ranks: " << rank_1 << "/" << size << "  &  ";
      std::cout << rank_2 << "/" << size << std::endl;
      double norm = full_pseudo_inverse.frobenius_norm();
      full_inverse.add(-1, full_pseudo_inverse);
      norm = full_inverse.frobenius_norm() / norm;
      AssertThrow(norm < tol, ExcInternalError());
    }
  pcout << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  const std::vector<unsigned int> sizes  = {{200, 400, 600}};
  const std::vector<unsigned int> blocks = {{32, 64}};

  const double tol = 1e-10;

  for (const auto &s : sizes)
    for (const auto &b : blocks)
      if (b <= s)
        test<double>(s, b, tol);
}
