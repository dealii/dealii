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

// test eigenpairs_symmetric_by_index(const std::pair<unsigned int,unsigned int>
// &, const bool) for some eigenvalues without eigenvectors


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

  pcout << size << ' ' << block_size << ' ' << grid->get_process_grid_rows()
        << ' ' << grid->get_process_grid_columns() << std::endl;

  const unsigned int n_eigenvalues     = size;
  const unsigned int max_n_eigenvalues = 5;

  // Create SPD matrices of requested size:
  FullMatrix<NumberType>  full_A(size);
  std::vector<NumberType> eigenvalues_Lapack(size);

  ScaLAPACKMatrix<NumberType> scalapack_syevx(size, grid, block_size);
  scalapack_syevx.set_property(LAPACKSupport::Property::symmetric);

  create_spd(full_A);
  scalapack_syevx = full_A;

  // Lapack as reference
  {
    std::vector<NumberType> lapack_A(size * size);
    for (unsigned int i = 0; i < size; ++i)
      for (unsigned int j = 0; j < size; ++j)
        lapack_A[i * size + j] = full_A(i, j);

    int info; // Variable containing information about the successful exit of
              // the lapack routine
    char jobz = 'V'; //'V': all eigenpairs of A are computed
    char uplo = 'U'; // storage format of the matrix A; not so important as
                     // matrix is symmetric
    int                     LDA = size; // leading dimension of the matrix A
    int                     lwork;      // length of vector/array work
    std::vector<NumberType> work(1);

    // by setting lwork to -1 a workspace query for work is done
    // as matrix is symmetric: LDA == size of matrix
    lwork = -1;
    syev(&jobz,
         &uplo,
         &LDA,
         &*lapack_A.begin(),
         &LDA,
         &*eigenvalues_Lapack.begin(),
         &*work.begin(),
         &lwork,
         &info);
    lwork = static_cast<int>(work[0]);
    work.resize(lwork);
    syev(&jobz,
         &uplo,
         &LDA,
         &*lapack_A.begin(),
         &LDA,
         &*eigenvalues_Lapack.begin(),
         &*work.begin(),
         &lwork,
         &info);
    AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("syev", info));
  }

  // the actual test:

  pcout << "comparing " << max_n_eigenvalues
        << " eigenvalues computed using LAPACK and ScaLAPACK pdsyevx:"
        << std::endl;
  const std::vector<NumberType> eigenvalues_psyevx =
    scalapack_syevx.eigenpairs_symmetric_by_index(
      std::make_pair(size - max_n_eigenvalues, size - 1), false);
  for (unsigned int i = eigenvalues_psyevx.size() - 1; i > 0; --i)
    {
      if (!(std::abs(eigenvalues_psyevx[i] -
                     eigenvalues_Lapack[size - eigenvalues_psyevx.size() + i]) /
              std::abs(
                eigenvalues_Lapack[size - eigenvalues_psyevx.size() + i]) <
            tol))
        {
          std::cout << "process #" << this_mpi_process
                    << ": eigenvalues do not fit: " << eigenvalues_psyevx[i]
                    << " <--> "
                    << eigenvalues_Lapack[size - eigenvalues_psyevx.size() + i]
                    << std::endl;
        }

      AssertThrow(
        std::abs(eigenvalues_psyevx[i] -
                 eigenvalues_Lapack[size - eigenvalues_psyevx.size() + i]) /
            std::abs(eigenvalues_Lapack[size - eigenvalues_psyevx.size() + i]) <
          tol,
        ExcInternalError());
    }
  pcout << "   with respect to the given tolerance the eigenvalues coincide"
        << std::endl;
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
