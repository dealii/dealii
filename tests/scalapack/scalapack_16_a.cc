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

// test pseudo_inverse(const NumberType ratio) for rectangular matrices

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
test(const unsigned int size_1,
     const unsigned int size_2,
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
      mpi_communicator, size_1, size_2, block_size, block_size);

  pcout << size_1 << "x" << size_2 << " " << block_size << " "
        << grid->get_process_grid_rows() << " "
        << grid->get_process_grid_columns() << std::endl;

  // Create randomly filled matrix of requested dimension size_1xsize_2:
  FullMatrix<NumberType> full_A(size_1, size_2);
  create_random(full_A);
  const NumberType ratio = 1e-8;

  ScaLAPACKMatrix<NumberType> pseudoinverse_A(size_1,
                                              size_2,
                                              grid,
                                              block_size,
                                              block_size,
                                              LAPACKSupport::Property::general);
  pseudoinverse_A         = full_A;
  const unsigned int rank = pseudoinverse_A.pseudoinverse(ratio);

  // The pseudoinverse_A has to fulfill: A * A+ * A = A

  ScaLAPACKMatrix<NumberType> scalapack_A(size_1,
                                          size_2,
                                          grid,
                                          block_size,
                                          block_size,
                                          LAPACKSupport::Property::general);
  scalapack_A = full_A;
  ScaLAPACKMatrix<NumberType> scalapack_tmp(size_2,
                                            size_2,
                                            grid,
                                            block_size,
                                            block_size,
                                            LAPACKSupport::Property::general);
  // compute tmp = A+ * A
  pseudoinverse_A.mmult(scalapack_tmp, scalapack_A);
  ScaLAPACKMatrix<NumberType> scalapack_result(
    size_1,
    size_2,
    grid,
    block_size,
    block_size,
    LAPACKSupport::Property::general);
  // compute result = A * tmp
  scalapack_A.mmult(scalapack_result, scalapack_tmp);

  // compute difference
  scalapack_result.add(-1, scalapack_A);
  const double norm = scalapack_result.frobenius_norm();

  if (this_mpi_process == 0)
    {
      std::cout << "rank=" << rank << "/" << std::min(size_1, size_2)
                << std::endl;
      AssertThrow(norm < tol, ExcInternalError());
    }
  pcout << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  const std::vector<unsigned int> sizes_1 = {{200, 400}};
  const std::vector<unsigned int> sizes_2 = {{250, 350}};
  const std::vector<unsigned int> blocks  = {{16, 32}};

  const double tol = 1e-10;

  for (const auto &s_1 : sizes_1)
    for (const auto &s_2 : sizes_2)
      for (const auto &b : blocks)
        if (b <= s_1 && b <= s_2)
          test<double>(s_1, s_2, b, tol);
}
