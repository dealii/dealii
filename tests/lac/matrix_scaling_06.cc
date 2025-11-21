// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test matrix scaling in distributed setting with trilinos and petsc matrices

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/matrix_scaling.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "../tests.h"

using namespace dealii;

template <typename MatrixType>
void
print_matrix(const MatrixType &matrix)
{
  const unsigned int this_process =
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_processes =
    Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  for (unsigned int proc = 0; proc < n_processes; ++proc)
    {
      if (proc == this_process)
        {
          matrix.print(deallog.get_file_stream());
        }
      MPI_Barrier(MPI_COMM_WORLD);
    }
}


template <typename MatrixType>
void
test(unsigned int scaling_steps)
{
  const unsigned int size = 10;

  bool converged;

  IndexSet locally_owned_dofs =
    Utilities::MPI::create_evenly_distributed_partitioning(MPI_COMM_WORLD,
                                                           size);
  IndexSet locally_relevant_dofs = locally_owned_dofs;

  DynamicSparsityPattern dsp(locally_relevant_dofs);

  for (auto row : locally_owned_dofs)
    {
      dsp.add(row, row);
    }
  dsp.add(0, size - 1);
  dsp.add(size - 1, 0);

  MatrixType matrix;
  matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, MPI_COMM_WORLD);

  for (auto row : locally_owned_dofs)
    {
      matrix.set(row, row, static_cast<double>(row + 1));
    }
  matrix.set(size - 1, size - 1, 1);
  matrix.set(0, size - 1, 1);
  matrix.set(size - 1, 0, 1);

  matrix.compress(VectorOperation::insert);

  print_matrix(matrix);

  MatrixScaling::AdditionalData additional_data;
  additional_data.l1linf_parameters.l1_norm_steps = scaling_steps;
  MatrixScaling scaler(additional_data);
  converged = scaler.find_scaling_and_scale_matrix(matrix);

  const Vector<double> &row_scaling    = scaler.get_row_scaling();
  const Vector<double> &column_scaling = scaler.get_column_scaling();

  deallog << "Reciprocal of scaling vectors squared " << std::endl;

  for (unsigned int i = 0; i < row_scaling.size(); i++)
    deallog << 1.0 / (row_scaling[i] * row_scaling[i]) << " ";
  deallog << std::endl;
  for (unsigned int i = 0; i < column_scaling.size(); i++)
    deallog << 1.0 / (column_scaling[i] * column_scaling[i]) << " ";
  deallog << std::endl;

  deallog << "Scaled matrix" << std::endl;
  print_matrix(matrix);

  deallog << "Converged? " << (converged ? "true" : "false") << std::endl;
}

template <typename MatrixType>
void
test_diagonal(unsigned int scaling_steps)
{
  const unsigned int size = 10;

  bool converged;

  IndexSet locally_owned_dofs =
    Utilities::MPI::create_evenly_distributed_partitioning(MPI_COMM_WORLD,
                                                           size);
  IndexSet locally_relevant_dofs = locally_owned_dofs;

  DynamicSparsityPattern dsp(locally_relevant_dofs);

  for (auto row : locally_owned_dofs)
    {
      dsp.add(row, row);
    }

  MatrixType matrix;
  matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, MPI_COMM_WORLD);

  for (auto row : locally_owned_dofs)
    {
      matrix.set(row, row, static_cast<double>(row + 1));
    }

  matrix.compress(VectorOperation::insert);

  print_matrix(matrix);

  MatrixScaling::AdditionalData additional_data;
  additional_data.l1linf_parameters.l1_norm_steps = scaling_steps;
  additional_data.algorithm =
    MatrixScaling::AdditionalData::ScalingAlgorithm::sinkhorn_knopp;
  additional_data.sinkhorn_knopp_parameters.norm_type =
    MatrixScaling::AdditionalData::SKParameters::NormType::l_infty;
  MatrixScaling scaler(additional_data);
  converged = scaler.find_scaling_and_scale_matrix(matrix);

  const Vector<double> &row_scaling    = scaler.get_row_scaling();
  const Vector<double> &column_scaling = scaler.get_column_scaling();

  deallog << "Reciprocal of scaling vectors squared " << std::endl;

  for (unsigned int i = 0; i < row_scaling.size(); i++)
    deallog << 1.0 / (row_scaling[i] * row_scaling[i]) << " ";
  deallog << std::endl;
  for (unsigned int i = 0; i < column_scaling.size(); i++)
    deallog << 1.0 / (column_scaling[i] * column_scaling[i]) << " ";
  deallog << std::endl;

  deallog << "Scaled matrix" << std::endl;
  print_matrix(matrix);

  deallog << "Converged? " << (converged ? "true" : "false") << std::endl;
}

template <typename MatrixType, typename VectorType>
void
test_scaling_linsys()
{
  const unsigned int size = 10;

  IndexSet locally_owned_dofs =
    Utilities::MPI::create_evenly_distributed_partitioning(MPI_COMM_WORLD,
                                                           size);

  DynamicSparsityPattern dsp(locally_owned_dofs);

  for (auto row : locally_owned_dofs)
    {
      dsp.add(row, row);
    }

  MatrixType matrix;
  matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, MPI_COMM_WORLD);

  for (auto row : locally_owned_dofs)
    {
      matrix.set(row, row, 1.0 + static_cast<double>(row));
    }
  matrix.compress(VectorOperation::insert);

  VectorType rhs, x, x_scaled;
  rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  x.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  x_scaled.reinit(locally_owned_dofs, MPI_COMM_WORLD);

  for (auto i : locally_owned_dofs)
    {
      rhs[i] = 1.0;
    }
  rhs.compress(VectorOperation::insert);

  if constexpr (std::is_same_v<MatrixType, TrilinosWrappers::SparseMatrix>)
    {
      TrilinosWrappers::SolverDirect solver;
      solver.initialize(matrix);
      solver.solve(x, rhs);
    }
  else if constexpr (std::is_same_v<MatrixType,
                                    PETScWrappers::MPI::SparseMatrix>)
    {
      SolverControl                   solver_control(1000, 1e-12);
      PETScWrappers::SolverCG         solver(solver_control);
      PETScWrappers::PreconditionNone preconditioner;
      preconditioner.initialize(matrix);
      solver.solve(matrix, x, rhs, preconditioner);
    }

  deallog << "Original system solution:" << std::endl;
  for (auto i : locally_owned_dofs)
    {
      deallog << x[i] << " ";
    }
  deallog << std::endl;

  MatrixType matrix_copy;
  matrix_copy.reinit(locally_owned_dofs,
                     locally_owned_dofs,
                     dsp,
                     MPI_COMM_WORLD);
  matrix_copy.copy_from(matrix);

  VectorType rhs_copy = rhs;

  MatrixScaling::AdditionalData control;
  MatrixScaling                 scaler(control);

  scaler.find_scaling_and_scale_linear_system(matrix_copy, rhs_copy);


  if constexpr (std::is_same_v<MatrixType, TrilinosWrappers::SparseMatrix>)
    {
      TrilinosWrappers::SolverDirect solver_scaled;
      solver_scaled.initialize(matrix_copy);
      solver_scaled.solve(x_scaled, rhs_copy);
    }
  else if constexpr (std::is_same_v<MatrixType,
                                    PETScWrappers::MPI::SparseMatrix>)
    {
      SolverControl                   solver_control_scaled(1000, 1e-12);
      PETScWrappers::SolverCG         solver_scaled(solver_control_scaled);
      PETScWrappers::PreconditionNone preconditioner_scaled;
      preconditioner_scaled.initialize(matrix_copy);
      solver_scaled.solve(matrix_copy,
                          x_scaled,
                          rhs_copy,
                          preconditioner_scaled);
    }

  scaler.scale_system_solution(x_scaled);


  deallog << "Scaled system solution (should match original):" << std::endl;
  for (auto i : locally_owned_dofs)
    {
      deallog << x_scaled[i] << " ";
    }
  deallog << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  MPILogInitAll log(true);

  unsigned int runs = 3;

  deallog << "Trilinos matrix" << std::endl;
  test<TrilinosWrappers::SparseMatrix>(runs);
  deallog << "Trilinos diagonal matrix " << std::endl;
  test_diagonal<TrilinosWrappers::SparseMatrix>(runs);
  deallog << "Trilinos linear system test" << std::endl;
  test_scaling_linsys<TrilinosWrappers::SparseMatrix,
                      TrilinosWrappers::MPI::Vector>();

  deallog << "PETSc matrix" << std::endl;
  test<PETScWrappers::MPI::SparseMatrix>(runs);
  deallog << "PETSc diagonal matrix " << std::endl;
  test_diagonal<PETScWrappers::MPI::SparseMatrix>(runs);
  deallog << "PETSc linear system test" << std::endl;
  test_scaling_linsys<PETScWrappers::MPI::SparseMatrix,
                      PETScWrappers::MPI::Vector>();
}
