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

// Verify that reusing an optimizer with changing constraint sets does not
// retain constraints from an earlier solve.

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/optimization/petsc_nonlinear_optimizer.templates.h>

#include <cmath>

#include "../tests.h"

using VectorType = PETScWrappers::MPI::Vector;
using MatrixType = PETScWrappers::MPI::SparseMatrix;
using Optimizer  = PETScWrappers::NonlinearOptimizer<VectorType, MatrixType>;

namespace
{
  void
  reinit_matrix(MatrixType &matrix)
  {
    Mat petsc_matrix;
    PetscCallAbort(
      MPI_COMM_SELF,
      MatCreateSeqAIJ(MPI_COMM_SELF, 1, 1, 1, nullptr, &petsc_matrix));
    static_cast<PETScWrappers::MatrixBase &>(matrix).reinit(petsc_matrix);
    PetscCallAbort(MPI_COMM_SELF, MatDestroy(&petsc_matrix));
  }



  void
  assert_solution(const VectorType &solution, const double expected_value)
  {
    AssertThrow(std::abs(solution(0) - expected_value) < 1e-5,
                ExcInternalError());
  }
} // namespace



int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  VectorType                            solution_vector(MPI_COMM_SELF, 1, 1);
  MatrixType                            hessian_matrix;
  Optimizer                             optimizer;
  std::optional<Optimizer::Constraints> equality_constraints;
  std::optional<Optimizer::Constraints> inequality_constraints;

  reinit_matrix(hessian_matrix);
  equality_constraints.emplace();
  equality_constraints->values.reinit(MPI_COMM_SELF, 1, 1);
  reinit_matrix(equality_constraints->jacobian);
  inequality_constraints.emplace();
  inequality_constraints->values.reinit(MPI_COMM_SELF, 1, 1);
  reinit_matrix(inequality_constraints->jacobian);

  optimizer.objective = [](const VectorType &x) -> PetscReal {
    return .5 * (x(0) - 2.) * (x(0) - 2.);
  };
  optimizer.gradient = [](const VectorType &x, VectorType &gradient_vector) {
    gradient_vector(0) = x(0) - 2.;
    gradient_vector.compress(VectorOperation::insert);
  };
  optimizer.hessian =
    [](const VectorType &, MatrixType &hessian_matrix, MatrixType &) {
      hessian_matrix = 0.;
      hessian_matrix.set(0, 0, 1.);
      hessian_matrix.compress(VectorOperation::insert);
    };
  optimizer.equality_constraints = [](const VectorType &x, VectorType &values) {
    values(0) = x(0) - 1.;
    values.compress(VectorOperation::insert);
  };
  optimizer.equality_jacobian =
    [](const VectorType &, MatrixType &matrix, MatrixType &) {
      matrix = 0.;
      matrix.set(0, 0, 1.);
      matrix.compress(VectorOperation::insert);
    };
  optimizer.inequality_constraints = [](const VectorType &x,
                                        VectorType       &values) {
    values(0) = x(0);
    values.compress(VectorOperation::insert);
  };
  optimizer.inequality_jacobian =
    [](const VectorType &, MatrixType &matrix, MatrixType &) {
      matrix = 0.;
      matrix.set(0, 0, 1.);
      matrix.compress(VectorOperation::insert);
    };

  solution_vector = 0.;
  optimizer.solve(solution_vector, hessian_matrix, equality_constraints);
  assert_solution(solution_vector, 1.);

  solution_vector = 0.;
  optimizer.solve(solution_vector,
                  hessian_matrix,
                  std::nullopt,
                  inequality_constraints);
  assert_solution(solution_vector, 2.);

  solution_vector = 0.;
  optimizer.solve(solution_vector, hessian_matrix);
  assert_solution(solution_vector, 2.);

  deallog << "OK" << std::endl;
}
