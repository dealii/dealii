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

// Solve min 1/2 ((x_0 - 1)^2 + (x_1 - 2)^2), whose minimizer is (1, 2),
// by configuring TAOALMM directly. The empty equality constraint and inactive
// variable bounds provide the constraint objects required by TAOALMM.

#include <deal.II/base/mpi.h>

#include <petsctao.h>

#include <cmath>

#include "../tests.h"

namespace
{
  PetscErrorCode
  objective(Tao, Vec x, PetscReal *value, void *)
  {
    const PetscScalar *entries;
    PetscFunctionBegin;
    PetscCall(VecGetArrayRead(x, &entries));
    *value =
      .5 *
      ((PetscRealPart(entries[0]) - 1.) * (PetscRealPart(entries[0]) - 1.) +
       (PetscRealPart(entries[1]) - 2.) * (PetscRealPart(entries[1]) - 2.));
    PetscCall(VecRestoreArrayRead(x, &entries));
    PetscFunctionReturn(PETSC_SUCCESS);
  }



  PetscErrorCode
  gradient(Tao, Vec x, Vec gradient, void *)
  {
    const PetscScalar *x_entries;
    PetscScalar       *gradient_entries;
    PetscFunctionBegin;
    PetscCall(VecGetArrayRead(x, &x_entries));
    PetscCall(VecGetArray(gradient, &gradient_entries));
    gradient_entries[0] = x_entries[0] - 1.;
    gradient_entries[1] = x_entries[1] - 2.;
    PetscCall(VecRestoreArray(gradient, &gradient_entries));
    PetscCall(VecRestoreArrayRead(x, &x_entries));
    PetscFunctionReturn(PETSC_SUCCESS);
  }



  PetscErrorCode
  hessian(Tao, Vec, Mat hessian, Mat, void *)
  {
    PetscFunctionBegin;
    PetscCall(MatZeroEntries(hessian));
    PetscCall(MatSetValue(hessian, 0, 0, 1., INSERT_VALUES));
    PetscCall(MatSetValue(hessian, 1, 1, 1., INSERT_VALUES));
    PetscCall(MatAssemblyBegin(hessian, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(hessian, MAT_FINAL_ASSEMBLY));
    PetscFunctionReturn(PETSC_SUCCESS);
  }



  PetscErrorCode
  equality_constraints(Tao, Vec, Vec, void *)
  {
    return PETSC_SUCCESS;
  }



  PetscErrorCode
  equality_jacobian(Tao, Vec, Mat jacobian, Mat, void *)
  {
    PetscFunctionBegin;
    PetscCall(MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
} // namespace



int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  Vec x;
  Vec gradient_vector;
  Vec equality_constraint_values;
  Vec lower_bounds;
  Vec upper_bounds;
  Mat hessian_matrix;
  Mat equality_jacobian_matrix;
  Tao tao;
  PetscCallAbort(MPI_COMM_SELF, VecCreateSeq(MPI_COMM_SELF, 2, &x));
  PetscCallAbort(MPI_COMM_SELF, VecSet(x, 0.));
  PetscCallAbort(MPI_COMM_SELF,
                 VecCreateSeq(MPI_COMM_SELF, 2, &gradient_vector));
  PetscCallAbort(MPI_COMM_SELF, VecSet(gradient_vector, 0.));
  PetscCallAbort(MPI_COMM_SELF,
                 VecCreateSeq(MPI_COMM_SELF, 0, &equality_constraint_values));
  PetscCallAbort(MPI_COMM_SELF, VecCreateSeq(MPI_COMM_SELF, 2, &lower_bounds));
  PetscCallAbort(MPI_COMM_SELF, VecCreateSeq(MPI_COMM_SELF, 2, &upper_bounds));
  PetscCallAbort(MPI_COMM_SELF, VecSet(lower_bounds, PETSC_NINFINITY));
  PetscCallAbort(MPI_COMM_SELF, VecSet(upper_bounds, PETSC_INFINITY));
  PetscCallAbort(
    MPI_COMM_SELF,
    MatCreateSeqAIJ(MPI_COMM_SELF, 2, 2, 2, nullptr, &hessian_matrix));
  PetscCallAbort(MPI_COMM_SELF,
                 MatSetValue(hessian_matrix, 0, 0, 0., INSERT_VALUES));
  PetscCallAbort(MPI_COMM_SELF,
                 MatSetValue(hessian_matrix, 1, 1, 0., INSERT_VALUES));
  PetscCallAbort(MPI_COMM_SELF,
                 MatAssemblyBegin(hessian_matrix, MAT_FINAL_ASSEMBLY));
  PetscCallAbort(MPI_COMM_SELF,
                 MatAssemblyEnd(hessian_matrix, MAT_FINAL_ASSEMBLY));
  PetscCallAbort(MPI_COMM_SELF,
                 MatCreateSeqAIJ(
                   MPI_COMM_SELF, 0, 2, 0, nullptr, &equality_jacobian_matrix));
  PetscCallAbort(MPI_COMM_SELF,
                 MatAssemblyBegin(equality_jacobian_matrix,
                                  MAT_FINAL_ASSEMBLY));
  PetscCallAbort(MPI_COMM_SELF,
                 MatAssemblyEnd(equality_jacobian_matrix, MAT_FINAL_ASSEMBLY));

  PetscCallAbort(MPI_COMM_SELF, TaoCreate(MPI_COMM_SELF, &tao));
  PetscCallAbort(MPI_COMM_SELF, TaoSetType(tao, TAOALMM));
#if DEAL_II_PETSC_VERSION_GTE(3, 17, 0)
  PetscCallAbort(MPI_COMM_SELF, TaoSetSolution(tao, x));
#else
  PetscCallAbort(MPI_COMM_SELF, TaoSetInitialVector(tao, x));
#endif
  PetscCallAbort(MPI_COMM_SELF,
                 TaoSetVariableBounds(tao, lower_bounds, upper_bounds));
#if DEAL_II_PETSC_VERSION_GTE(3, 17, 0)
  PetscCallAbort(MPI_COMM_SELF, TaoSetObjective(tao, objective, nullptr));
  PetscCallAbort(MPI_COMM_SELF,
                 TaoSetGradient(tao, gradient_vector, gradient, nullptr));
#else
  PetscCallAbort(MPI_COMM_SELF,
                 TaoSetObjectiveRoutine(tao, objective, nullptr));
  PetscCallAbort(MPI_COMM_SELF, TaoSetGradientRoutine(tao, gradient, nullptr));
#endif
  PetscCallAbort(MPI_COMM_SELF, TaoSetFromOptions(tao));
  PetscCallAbort(MPI_COMM_SELF,
                 TaoSetEqualityConstraintsRoutine(tao,
                                                  equality_constraint_values,
                                                  equality_constraints,
                                                  nullptr));
  PetscCallAbort(MPI_COMM_SELF,
                 TaoSetJacobianEqualityRoutine(tao,
                                               equality_jacobian_matrix,
                                               equality_jacobian_matrix,
                                               equality_jacobian,
                                               nullptr));

#if DEAL_II_PETSC_VERSION_GTE(3, 17, 0)
  PetscCallAbort(
    MPI_COMM_SELF,
    TaoSetHessian(tao, hessian_matrix, hessian_matrix, hessian, nullptr));
#else
  PetscCallAbort(MPI_COMM_SELF,
                 TaoSetHessianRoutine(
                   tao, hessian_matrix, hessian_matrix, hessian, nullptr));
#endif
  PetscCallAbort(MPI_COMM_SELF, TaoSolve(tao));
  PetscInt iterations;
  PetscCallAbort(MPI_COMM_SELF, TaoGetIterationNumber(tao, &iterations));
  deallog << "The optimizer found the solution in " << iterations
          << " iterations." << std::endl;

  const PetscScalar *solution;
  PetscCallAbort(MPI_COMM_SELF, VecGetArrayRead(x, &solution));
  AssertThrow(std::abs(PetscRealPart(solution[0]) - 1.) < 1e-6 &&
                std::abs(PetscRealPart(solution[1]) - 2.) < 1e-6,
              ExcInternalError());
  PetscCallAbort(MPI_COMM_SELF, VecRestoreArrayRead(x, &solution));

  PetscCallAbort(MPI_COMM_SELF, TaoDestroy(&tao));
  PetscCallAbort(MPI_COMM_SELF, MatDestroy(&equality_jacobian_matrix));
  PetscCallAbort(MPI_COMM_SELF, MatDestroy(&hessian_matrix));
  PetscCallAbort(MPI_COMM_SELF, VecDestroy(&upper_bounds));
  PetscCallAbort(MPI_COMM_SELF, VecDestroy(&lower_bounds));
  PetscCallAbort(MPI_COMM_SELF, VecDestroy(&equality_constraint_values));
  PetscCallAbort(MPI_COMM_SELF, VecDestroy(&gradient_vector));
  PetscCallAbort(MPI_COMM_SELF, VecDestroy(&x));
  deallog << "OK" << std::endl;
}
