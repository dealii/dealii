// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

// test norms of ScaLAPACK vs FullMatrix

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/multithread_info.h>

#include <deal.II/lac/vector.h>

#include <deal.II/lac/scalapack.h>

#include <fstream>
#include <iostream>

template <typename NumberType>
void test(const unsigned int size, const unsigned int block_size, const NumberType tol)
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);
  const unsigned int n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator));
  const unsigned int this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator));

  ConditionalOStream pcout (std::cout, (this_mpi_process ==0));

  // Create SPD matrices of requested size:
  FullMatrix<NumberType> full_A(size), full_B(size);

  std::shared_ptr<Utilities::MPI::ProcessGrid> grid = std::make_shared<Utilities::MPI::ProcessGrid>(mpi_communicator,size,size,block_size,block_size);
  ScaLAPACKMatrix<NumberType> scalapack_A (size, grid, block_size,
                                           LAPACKSupport::Property::symmetric);
  ScaLAPACKMatrix<NumberType> scalapack_B (size, size, grid, block_size, block_size,
                                           LAPACKSupport::Property::general);

  pcout << size << " " << block_size << " " << grid->get_process_grid_rows() << " " << grid->get_process_grid_columns() << std::endl;

  create_spd (full_A);
  create_random (full_B);

  const NumberType l1_A = full_A.l1_norm();
  const NumberType linfty_A = full_A.linfty_norm();
  const NumberType frobenius_A = full_A.frobenius_norm();

  const NumberType l1_B = full_B.l1_norm();
  const NumberType linfty_B = full_B.linfty_norm();
  const NumberType frobenius_B = full_B.frobenius_norm();

  scalapack_A = full_A;
  scalapack_B = full_B;

  // local result on this core:
  const NumberType s_l1_A = scalapack_A.l1_norm();
  const NumberType s_linfty_A = scalapack_A.linfty_norm();
  const NumberType s_frobenius_A = scalapack_A.frobenius_norm();

  const NumberType s_l1_B = scalapack_B.l1_norm();
  const NumberType s_linfty_B = scalapack_B.linfty_norm();
  const NumberType s_frobenius_B = scalapack_B.frobenius_norm();

  // make sure we have the same result on all cores, do average:
  const NumberType as_l1_A = dealii::Utilities::MPI::sum(s_l1_A, mpi_communicator)/n_mpi_processes;
  const NumberType as_linfty_A = dealii::Utilities::MPI::sum(s_linfty_A, mpi_communicator)/n_mpi_processes;
  const NumberType as_frobenius_A = dealii::Utilities::MPI::sum(s_frobenius_A, mpi_communicator)/n_mpi_processes;

  const NumberType as_l1_B = dealii::Utilities::MPI::sum(s_l1_B, mpi_communicator)/n_mpi_processes;
  const NumberType as_linfty_B = dealii::Utilities::MPI::sum(s_linfty_B, mpi_communicator)/n_mpi_processes;
  const NumberType as_frobenius_B = dealii::Utilities::MPI::sum(s_frobenius_B, mpi_communicator)/n_mpi_processes;

  pcout << l1_A << " " << s_l1_A << " " << as_l1_A << std::endl
        << linfty_A << " " << s_linfty_A << " " << as_linfty_A << std::endl
        << frobenius_A << " " << s_frobenius_A << " " << as_frobenius_A << std::endl;

  pcout << l1_B << " " << s_l1_B << " " << as_l1_B << std::endl
        << linfty_B << " " << s_linfty_B << " " << as_linfty_B << std::endl
        << frobenius_B << " " << s_frobenius_B << " " << as_frobenius_B << std::endl;

  AssertThrow (std::abs(l1_A -as_l1_A) < tol*std::abs(l1_A), dealii::ExcInternalError());
  AssertThrow (std::abs(linfty_A -as_linfty_A) < tol*std::abs(linfty_A), dealii::ExcInternalError());
  AssertThrow (std::abs(frobenius_A -as_frobenius_A) < tol*std::abs(frobenius_A), dealii::ExcInternalError());

  AssertThrow (std::abs(l1_B -as_l1_B) < tol*std::abs(l1_B), dealii::ExcInternalError());
  AssertThrow (std::abs(linfty_B -as_linfty_B) < tol*std::abs(linfty_B), dealii::ExcInternalError());
  AssertThrow (std::abs(frobenius_B -as_frobenius_B) < tol*std::abs(frobenius_B), dealii::ExcInternalError());
}



int main (int argc,char **argv)
{

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);

  const std::vector<unsigned int> sizes = {{32,64,120,320,640}};
  const std::vector<unsigned int> blocks = {{32,64}};

  const double tol_double = 1e-10;
  const float tol_float = 1e-5;

  /*for (const auto &s : sizes)
    for (const auto &b : blocks)
      if (b <= s)
        test<float>(s,b,tol_float);*/

  for (const auto &s : sizes)
    for (const auto &b : blocks)
      if (b <= s)
        test<double>(s,b,tol_double);
}
