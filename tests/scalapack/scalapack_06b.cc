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

// test eigenpairs_symmetric_by_index(const std::pair<unsigned int,unsigned int>
// &, const bool) for some eigenvalues with eigenvectors


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

  const unsigned int n_eigenvalues     = size;
  const unsigned int max_n_eigenvalues = 5;

  // Create SPD matrices of requested size:
  FullMatrix<NumberType>          full_A(size);
  std::vector<NumberType>         eigenvalues_Lapack(size);
  std::vector<Vector<NumberType>> s_eigenvectors_(max_n_eigenvalues,
                                                  Vector<NumberType>(size));
  std::vector<Vector<NumberType>> p_eigenvectors_(max_n_eigenvalues,
                                                  Vector<NumberType>(size));
  FullMatrix<NumberType>          p_eigenvectors(size, size);

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
    for (unsigned int i = 0; i < max_n_eigenvalues; ++i)
      for (unsigned int j = 0; j < size; ++j)
        s_eigenvectors_[i][j] = lapack_A[(size - 1 - i) * size + j];
  }

  // the actual test:

  pcout
    << "comparing " << max_n_eigenvalues
    << " eigenvalues and eigenvectors computed using LAPACK and ScaLAPACK pdsyevx:"
    << std::endl;
  const std::vector<NumberType> eigenvalues_psyevx =
    scalapack_syevx.eigenpairs_symmetric_by_index(
      std::make_pair(size - max_n_eigenvalues, size - 1), true);
  scalapack_syevx.copy_to(p_eigenvectors);
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

  // FIXME: run-time error on macOS if code between "pcout << comparing" and
  // this line is executed. Happens at the end of the program run while freeing
  // the MPI communicator.

  for (unsigned int i = 0; i < max_n_eigenvalues; ++i)
    for (unsigned int j = 0; j < size; ++j)
      p_eigenvectors_[i][j] = p_eigenvectors(j, max_n_eigenvalues - 1 - i);

  // product of eigenvectors computed using Lapack and ScaLapack has to be
  // either 1 or -1
  for (unsigned int i = 0; i < max_n_eigenvalues; ++i)
    {
      const NumberType product = p_eigenvectors_[i] * s_eigenvectors_[i];
      if (!(std::abs(std::abs(product) - 1) < tol * 10))
        std::cout << "process #" << this_mpi_process
                  << ": eigenvectors do not coincide: abs(" << product
                  << ") != 1" << std::endl;

      // the requirement for alignment of the eigenvectors has to be released
      // (primarily for floats)
      AssertThrow(std::abs(std::abs(product) - 1) < tol * 10,
                  ExcInternalError());
    }
  pcout
    << "   with respect to the given tolerance also the eigenvectors coincide"
    << std::endl
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
