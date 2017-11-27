// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include "../tests.h"
#include "../lapack/create_matrix.h"

// test eigenvalues_symmetric()

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/multithread_info.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include <deal.II/lac/scalapack.h>

#include <fstream>
#include <iostream>
#include <algorithm>

extern "C"  //Some Lapack routines
{
  void dsyev_(char *jobz, char *uplo, int *n, double *A, int *lda, double *w, double *work, int *lwork, int *info);
  void ssyev_(char *jobz, char *uplo, int *n, float *A, int *lda, float *w, float *work, int *lwork, int *info);
}

template <typename NumberType>
void test(const unsigned int size, const unsigned int block_size)
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);
  const unsigned int n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator));
  const unsigned int this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator));

  ConditionalOStream pcout (std::cout, (this_mpi_process ==0));

  // Create SPD matrices of requested size:
  FullMatrix<NumberType> full_A(size);
  std::vector<NumberType> lapack_A(size*size);

  std::shared_ptr<ProcessGrid> grid = std::make_shared<ProcessGrid>(mpi_communicator,size,size,block_size,block_size);
  ScaLAPACKMatrix<NumberType> scalapack_A (size, grid, block_size,
                                           LAPACKSupport::Property::symmetric);

  pcout << size << " " << block_size << " " << grid->get_process_grid_rows() << " " << grid->get_process_grid_columns() << std::endl;
  create_spd (full_A);

  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      lapack_A[i*size+j] = full_A(i,j);

  std::vector<NumberType> eigenvalues_ScaLapack, eigenvalues_Lapack(size);
  //Lapack as reference
  {
    int info; //Variable containing information about the successfull exit of the lapack routine
    char jobz = 'N';  //'N': the eigenvalues_Lapack of A are computed
    char uplo = 'U';  //storage format of the matrix A; not so important as matrix is symmetric
    int LDA = size;   //leading dimension of the matrix A
    int lwork;      //length of vector/array work
    std::vector<double> work (1);

    //by setting lwork to -1 a workspace query for work is done
    //as matrix is symmetric: LDA == size of matrix
    lwork = -1;
    dsyev_(&jobz, &uplo, &LDA, & *lapack_A.begin(), &LDA, & *eigenvalues_Lapack.begin(), & *work.begin(), &lwork, &info);
    lwork=work[0];
    work.resize (lwork);
    dsyev_(&jobz, &uplo, &LDA, & *lapack_A.begin(), &LDA, & *eigenvalues_Lapack.begin(), & *work.begin(), &lwork, &info);

    AssertThrow (info==0, LAPACKSupport::ExcErrorCode("syev", info));
  }
  // Scalapack:
  scalapack_A = full_A;
  eigenvalues_ScaLapack = scalapack_A.eigenvalues_symmetric();
  unsigned int n_eigenvalues = eigenvalues_ScaLapack.size(), max_n_eigenvalues=5;

  pcout << "First " << max_n_eigenvalues << " ScaLapack eigenvalues" << std::endl;

  for (unsigned int i=0; i<max_n_eigenvalues; ++i)
    pcout << eigenvalues_ScaLapack[i] << "   ";

  pcout << std::endl << "First " << max_n_eigenvalues << " Lapack eigenvalues" << std::endl;

  for (unsigned int i=0; i<max_n_eigenvalues; ++i)
    pcout << eigenvalues_Lapack[i] << "   ";

  for (unsigned int i=0; i<max_n_eigenvalues; ++i)
    AssertThrow ( std::abs(eigenvalues_ScaLapack[i]-eigenvalues_Lapack[i]) < std::abs(eigenvalues_Lapack[i])*1e-10, dealii::ExcInternalError());

  pcout << std::endl << std::endl;
}



int main (int argc,char **argv)
{

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);

  const std::vector<unsigned int> sizes = {{32,64,120,320,640}};
  const std::vector<unsigned int> blocks = {{32,64}};

  for (const auto &s : sizes)
    for (const auto &b : blocks)
      if (b <= s)
        test<double>(s,b);

}
