// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Solve min sum_{i=0}^{19} x_i, subject to x_i>=0 and x_{i+1}-x_i>=1 for
// i=0,...,18. The unique minimizer is x_i=i. The initial guess is
// x_i=i+10+U(0,100).

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/optimization/petsc_nonlinear_optimizer.templates.h>

#include <cmath>
#include <random>

#include "../tests.h"

using VectorType = PETScWrappers::MPI::Vector;
using MatrixType = PETScWrappers::MPI::SparseMatrix;
using Optimizer  = PETScWrappers::NonlinearOptimizer<VectorType, MatrixType>;

namespace
{
  void
  reinit_matrix(MatrixType    &matrix,
                const PetscInt rows,
                const PetscInt columns,
                const PetscInt entries_per_row)
  {
    Mat petsc_matrix;
    PetscCallAbort(
      MPI_COMM_SELF,
      MatCreateSeqAIJ(
        MPI_COMM_SELF, rows, columns, entries_per_row, nullptr, &petsc_matrix));
    static_cast<PETScWrappers::MatrixBase &>(matrix).reinit(petsc_matrix);
    PetscCallAbort(MPI_COMM_SELF, MatDestroy(&petsc_matrix));
  }
} // namespace



int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  constexpr PetscInt n_variables    = 20;
  constexpr PetscInt n_inequalities = n_variables - 1;

  VectorType solution_vector(MPI_COMM_SELF, n_variables, n_variables);
  MatrixType hessian_matrix;
  reinit_matrix(hessian_matrix, n_variables, n_variables, 1);

  std::optional<Optimizer::Constraints> inequalities;
  inequalities.emplace();
  inequalities->values.reinit(MPI_COMM_SELF, n_inequalities, n_inequalities);
  reinit_matrix(inequalities->jacobian, n_inequalities, n_variables, 2);

  std::optional<Optimizer::Bounds> bounds;
  bounds.emplace();
  bounds->lower.reinit(MPI_COMM_SELF, n_variables, n_variables);
  bounds->upper.reinit(MPI_COMM_SELF, n_variables, n_variables);
  for (PetscInt index = 0; index < n_variables; ++index)
    {
      std::mt19937                           generator(index);
      std::uniform_real_distribution<double> distribution(0., 100.);
      solution_vector(index) = index + 10. + distribution(generator);
      bounds->lower(index)   = 0.;
      bounds->upper(index)   = PETSC_INFINITY;
    }
  solution_vector.compress(VectorOperation::insert);
  bounds->lower.compress(VectorOperation::insert);
  bounds->upper.compress(VectorOperation::insert);

  Optimizer::AdditionalData additional_data;
  additional_data.almm_subsolver =
    Optimizer::AdditionalData::ALMMSubsolver::bncg;
  Optimizer optimizer(additional_data);
  optimizer.objective = [](const VectorType &x) -> PetscReal {
    PetscReal result = 0.;
    for (PetscInt index = 0; index < n_variables; ++index)
      result += x(index);
    return result;
  };
  optimizer.gradient = [](const VectorType &, VectorType &gradient_vector) {
    gradient_vector = 1.;
  };
  optimizer.hessian =
    [](const VectorType &, MatrixType &hessian, MatrixType &) { hessian = 0.; };
  optimizer.inequality_constraints = [](const VectorType &x,
                                        VectorType       &values) {
    for (PetscInt index = 0; index < n_inequalities; ++index)
      values(index) = x(index + 1) - x(index) - 1.;
    values.compress(VectorOperation::insert);
  };
  optimizer.inequality_jacobian =
    [](const VectorType &, MatrixType &jacobian, MatrixType &) {
      jacobian = 0.;
      for (PetscInt index = 0; index < n_inequalities; ++index)
        {
          jacobian.set(index, index, -1.);
          jacobian.set(index, index + 1, 1.);
        }
      jacobian.compress(VectorOperation::insert);
    };

  const auto iterations = optimizer.solve(
    solution_vector, hessian_matrix, std::nullopt, inequalities, bounds);
  deallog << "The optimizer found the solution in " << iterations
          << " iterations." << std::endl;

  for (PetscInt index = 0; index < n_variables; ++index)
    AssertThrow(std::abs(solution_vector(index) - index) < 1e-6,
                ExcInternalError());
  deallog << "OK" << std::endl;
}
