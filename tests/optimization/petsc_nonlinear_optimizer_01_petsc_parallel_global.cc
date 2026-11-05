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

// On four MPI processes, directly configure TAOALMM to solve
// min 1/2 sum_{i=0}^7 (x_i-(i+1))^2 + 1/2 (sum_{i=0}^7 x_i-36)^2.
// The global coupling term leaves x_i=i+1 as the minimizer.

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
    PetscReal          local_value;
    PetscFunctionBegin;
    PetscCall(VecGetArrayRead(x, &entries));
    const PetscInt first =
      2 * Utilities::MPI::this_mpi_process(PetscObjectComm((PetscObject)x));
    local_value         = .5 * ((PetscRealPart(entries[0]) - (first + 1.)) *
                          (PetscRealPart(entries[0]) - (first + 1.)) +
                        (PetscRealPart(entries[1]) - (first + 2.)) *
                          (PetscRealPart(entries[1]) - (first + 2.)));
    PetscReal local_sum = PetscRealPart(entries[0]) + PetscRealPart(entries[1]);
    PetscReal sum;
    PetscCall(VecRestoreArrayRead(x, &entries));
    PetscCallMPI(MPIU_Allreduce(&local_sum,
                                &sum,
                                1,
                                MPIU_REAL,
                                MPIU_SUM,
                                PetscObjectComm((PetscObject)x)));
    PetscReal global_value;
    PetscCallMPI(MPIU_Allreduce(&local_value,
                                &global_value,
                                1,
                                MPIU_REAL,
                                MPIU_SUM,
                                PetscObjectComm((PetscObject)x)));
    *value = global_value + .5 * (sum - 36.) * (sum - 36.);
    PetscFunctionReturn(PETSC_SUCCESS);
  }



  PetscErrorCode
  gradient(Tao, Vec x, Vec gradient, void *)
  {
    const PetscScalar *x_entries;
    PetscScalar       *gradient_entries;
    PetscFunctionBegin;
    const MPI_Comm comm  = PetscObjectComm((PetscObject)x);
    const PetscInt first = 2 * Utilities::MPI::this_mpi_process(comm);
    PetscCall(VecGetArrayRead(x, &x_entries));
    PetscReal local_sum =
      PetscRealPart(x_entries[0]) + PetscRealPart(x_entries[1]);
    PetscReal sum;
    PetscCallMPI(
      MPIU_Allreduce(&local_sum, &sum, 1, MPIU_REAL, MPIU_SUM, comm));
    PetscCall(VecGetArray(gradient, &gradient_entries));
    gradient_entries[0] = x_entries[0] - (first + 1.) + sum - 36.;
    gradient_entries[1] = x_entries[1] - (first + 2.) + sum - 36.;
    PetscCall(VecRestoreArray(gradient, &gradient_entries));
    PetscCall(VecRestoreArrayRead(x, &x_entries));
    PetscFunctionReturn(PETSC_SUCCESS);
  }



  PetscErrorCode
  hessian(Tao, Vec, Mat hessian, Mat, void *)
  {
    PetscInt begin, end;
    PetscFunctionBegin;
    PetscCall(MatZeroEntries(hessian));
    PetscCall(MatGetOwnershipRange(hessian, &begin, &end));
    for (PetscInt row = begin; row < end; ++row)
      for (PetscInt column = 0; column < 8; ++column)
        PetscCall(MatSetValue(
          hessian, row, column, row == column ? 2. : 1., INSERT_VALUES));
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

  const MPI_Comm comm = MPI_COMM_WORLD;
  Vec            x, gradient_vector, equality_constraint_values, lower_bounds,
    upper_bounds;
  Mat hessian_matrix, equality_jacobian_matrix;
  Tao tao;
  PetscCallAbort(comm, VecCreateMPI(comm, 2, 8, &x));
  PetscCallAbort(comm, VecSet(x, 0.));
  PetscCallAbort(comm, VecCreateMPI(comm, 2, 8, &gradient_vector));
  PetscCallAbort(comm, VecCreateMPI(comm, 0, 0, &equality_constraint_values));
  PetscCallAbort(comm, VecCreateMPI(comm, 2, 8, &lower_bounds));
  PetscCallAbort(comm, VecCreateMPI(comm, 2, 8, &upper_bounds));
  PetscCallAbort(comm, VecSet(lower_bounds, PETSC_NINFINITY));
  PetscCallAbort(comm, VecSet(upper_bounds, PETSC_INFINITY));
  PetscCallAbort(
    comm,
    MatCreateAIJ(comm, 2, 2, 8, 8, 2, nullptr, 6, nullptr, &hessian_matrix));
  PetscCallAbort(comm, MatAssemblyBegin(hessian_matrix, MAT_FINAL_ASSEMBLY));
  PetscCallAbort(comm, MatAssemblyEnd(hessian_matrix, MAT_FINAL_ASSEMBLY));
  PetscCallAbort(
    comm,
    MatCreateAIJ(
      comm, 0, 2, 0, 8, 0, nullptr, 0, nullptr, &equality_jacobian_matrix));
  PetscCallAbort(comm,
                 MatAssemblyBegin(equality_jacobian_matrix,
                                  MAT_FINAL_ASSEMBLY));
  PetscCallAbort(comm,
                 MatAssemblyEnd(equality_jacobian_matrix, MAT_FINAL_ASSEMBLY));

  PetscCallAbort(comm, TaoCreate(comm, &tao));
  PetscCallAbort(comm, TaoSetType(tao, TAOALMM));
#if DEAL_II_PETSC_VERSION_GTE(3, 17, 0)
  PetscCallAbort(comm, TaoSetSolution(tao, x));
#else
  PetscCallAbort(comm, TaoSetInitialVector(tao, x));
#endif
  PetscCallAbort(comm, TaoSetVariableBounds(tao, lower_bounds, upper_bounds));
#if DEAL_II_PETSC_VERSION_GTE(3, 17, 0)
  PetscCallAbort(comm, TaoSetObjective(tao, objective, nullptr));
  PetscCallAbort(comm, TaoSetGradient(tao, gradient_vector, gradient, nullptr));
#else
  PetscCallAbort(comm, TaoSetObjectiveRoutine(tao, objective, nullptr));
  PetscCallAbort(comm, TaoSetGradientRoutine(tao, gradient, nullptr));
#endif
  PetscCallAbort(comm,
                 TaoSetEqualityConstraintsRoutine(tao,
                                                  equality_constraint_values,
                                                  equality_constraints,
                                                  nullptr));
  PetscCallAbort(comm,
                 TaoSetJacobianEqualityRoutine(tao,
                                               equality_jacobian_matrix,
                                               equality_jacobian_matrix,
                                               equality_jacobian,
                                               nullptr));
#if DEAL_II_PETSC_VERSION_GTE(3, 17, 0)
  PetscCallAbort(
    comm, TaoSetHessian(tao, hessian_matrix, hessian_matrix, hessian, nullptr));
#else
  PetscCallAbort(comm,
                 TaoSetHessianRoutine(
                   tao, hessian_matrix, hessian_matrix, hessian, nullptr));
#endif
  PetscCallAbort(comm, TaoSetFromOptions(tao));
  PetscCallAbort(comm, TaoSolve(tao));
  PetscInt iterations;
  PetscCallAbort(comm, TaoGetIterationNumber(tao, &iterations));
  deallog << "The optimizer found the solution in " << iterations
          << " iterations." << std::endl;

  const PetscScalar *solution;
  PetscCallAbort(comm, VecGetArrayRead(x, &solution));
  const PetscInt first = 2 * Utilities::MPI::this_mpi_process(comm);
  AssertThrow(std::abs(PetscRealPart(solution[0]) - (first + 1.)) < 1e-4 &&
                std::abs(PetscRealPart(solution[1]) - (first + 2.)) < 1e-4,
              ExcInternalError());
  PetscCallAbort(comm, VecRestoreArrayRead(x, &solution));

  PetscCallAbort(comm, TaoDestroy(&tao));
  PetscCallAbort(comm, MatDestroy(&equality_jacobian_matrix));
  PetscCallAbort(comm, MatDestroy(&hessian_matrix));
  PetscCallAbort(comm, VecDestroy(&upper_bounds));
  PetscCallAbort(comm, VecDestroy(&lower_bounds));
  PetscCallAbort(comm, VecDestroy(&equality_constraint_values));
  PetscCallAbort(comm, VecDestroy(&gradient_vector));
  PetscCallAbort(comm, VecDestroy(&x));
  deallog << "OK" << std::endl;
}
