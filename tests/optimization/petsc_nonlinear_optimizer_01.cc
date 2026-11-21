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

// Solve min 1/2 ((x_0 - 1)^2 + (x_1 - 2)^2), whose unconstrained minimizer
// is (1, 2), using PETScWrappers::NonlinearOptimizer.

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/optimization/petsc_linear_optimizer.h>
#include <deal.II/optimization/petsc_nonlinear_optimizer.templates.h>

#include <cmath>

#include "../tests.h"

using VectorType = PETScWrappers::MPI::Vector;
using MatrixType = PETScWrappers::MPI::SparseMatrix;
using Optimizer  = PETScWrappers::NonlinearOptimizer<VectorType, MatrixType>;

namespace
{
  void
  reinit_matrix(MatrixType &matrix, const PetscInt rows, const PetscInt columns)
  {
    Mat petsc_matrix;
    PetscCallAbort(
      MPI_COMM_SELF,
      MatCreateSeqAIJ(
        MPI_COMM_SELF, rows, columns, columns, nullptr, &petsc_matrix));
    static_cast<PETScWrappers::MatrixBase &>(matrix).reinit(petsc_matrix);
    PetscCallAbort(MPI_COMM_SELF, MatDestroy(&petsc_matrix));
  }



  void
  set_identity(MatrixType &matrix)
  {
    matrix = 0.;
    matrix.set(0, 0, 1.);
    matrix.set(1, 1, 1.);
    matrix.compress(VectorOperation::insert);
  }



  void
  set_quadratic_callbacks(Optimizer   &optimizer,
                          const double target_0,
                          const double target_1)
  {
    optimizer.objective = [=](const VectorType &x) -> PetscReal {
      return .5 * ((x(0) - target_0) * (x(0) - target_0) +
                   (x(1) - target_1) * (x(1) - target_1));
    };
    optimizer.gradient = [=](const VectorType &x, VectorType &gradient_vector) {
      gradient_vector(0) = x(0) - target_0;
      gradient_vector(1) = x(1) - target_1;
      gradient_vector.compress(VectorOperation::insert);
    };
    optimizer.hessian = [](const VectorType &,
                           MatrixType &hessian_matrix,
                           MatrixType &) { set_identity(hessian_matrix); };
  }



  void
  assert_solution(const VectorType &solution,
                  const double      expected_0,
                  const double      expected_1)
  {
    AssertThrow(std::abs(solution(0) - expected_0) < 1e-6 &&
                  std::abs(solution(1) - expected_1) < 1e-6,
                ExcInternalError());
  }



  void
  test_unconstrained()
  {
    VectorType solution_vector(MPI_COMM_SELF, 2, 2);
    MatrixType hessian_matrix;
    Optimizer  optimizer;

    reinit_matrix(hessian_matrix, 2, 2);
    solution_vector = 0.;
    set_quadratic_callbacks(optimizer, 1., 2.);
    const auto iterations = optimizer.solve(solution_vector, hessian_matrix);
    deallog << "The optimizer found the solution in " << iterations
            << " iterations." << std::endl;
    assert_solution(solution_vector, 1., 2.);
  }



} // namespace



int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  test_unconstrained();

  deallog << "OK" << std::endl;
}
