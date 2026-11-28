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

// Solve min x_0+2x_1, subject to x_0+x_1=3 and x_0<=2. The minimizer is
// (2,1).

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/optimization/petsc_linear_optimizer.h>

#include <cmath>

#include "../tests.h"

using VectorType = PETScWrappers::MPI::Vector;
using MatrixType = PETScWrappers::MPI::SparseMatrix;

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
} // namespace



int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  VectorType solution_vector(MPI_COMM_SELF, 2, 2);
  VectorType objective(MPI_COMM_SELF, 2, 2);
  VectorType equality_right_hand_side(MPI_COMM_SELF, 1, 1);
  VectorType inequality_right_hand_side(MPI_COMM_SELF, 1, 1);
  MatrixType equality_matrix;
  MatrixType inequality_matrix;
  reinit_matrix(equality_matrix, 1, 2);
  reinit_matrix(inequality_matrix, 1, 2);

  solution_vector(0) = 1.;
  solution_vector(1) = 2.;
  solution_vector.compress(VectorOperation::insert);
  objective(0) = 1.;
  objective(1) = 2.;
  objective.compress(VectorOperation::insert);
  equality_matrix.set(0, 0, 1.);
  equality_matrix.set(0, 1, 1.);
  equality_matrix.compress(VectorOperation::insert);
  equality_right_hand_side(0) = 3.;
  equality_right_hand_side.compress(VectorOperation::insert);
  inequality_matrix.set(0, 0, -1.);
  inequality_matrix.compress(VectorOperation::insert);
  inequality_right_hand_side(0) = -2.;
  inequality_right_hand_side.compress(VectorOperation::insert);

  PETScWrappers::LinearOptimizer<VectorType, MatrixType> optimizer;
  const auto iterations = optimizer.solve(solution_vector,
                                          objective,
                                          &equality_matrix,
                                          &equality_right_hand_side,
                                          &inequality_matrix,
                                          &inequality_right_hand_side);
  deallog << "The optimizer found the solution in " << iterations
          << " iterations." << std::endl;

  AssertThrow(std::abs(solution_vector(0) - 2.) < 1e-6 &&
                std::abs(solution_vector(1) - 1.) < 1e-6,
              ExcInternalError());
  deallog << "OK" << std::endl;
}
