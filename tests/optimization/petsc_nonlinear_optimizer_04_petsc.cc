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

// Solve min sum_{i=0}^{n-1} x_i, subject to x_i>=0 and x_{j+1}-x_j>=1 for
// j=0,...,n-2, by configuring TAOALMM directly. The unique minimizer is x_i=i.
// The initial guess is x_i=i+10+U(0,100).

#include <petsctao.h>

#include <cmath>
#include <iostream>
#include <random>

#if !PETSC_VERSION_GE(3, 19, 0)
#  define PETSC_SUCCESS 0
#endif

#if !PETSC_VERSION_GE(3, 17, 0)
#  define PetscCall(expression)                               \
    do                                                        \
      {                                                       \
        const PetscErrorCode petsc_error_code = (expression); \
        if (petsc_error_code != 0)                            \
          return petsc_error_code;                            \
      }                                                       \
    while (false)
#  define PetscCallAbort(comm, expression)                    \
    do                                                        \
      {                                                       \
        const PetscErrorCode petsc_error_code = (expression); \
        if (petsc_error_code != 0)                            \
          {                                                   \
            MPI_Abort(comm, petsc_error_code);                \
            return petsc_error_code;                          \
          }                                                   \
      }                                                       \
    while (false)
#endif

namespace
{
  constexpr PetscInt n_variables    = 100;
  constexpr PetscInt n_inequalities = n_variables - 1;

  PetscErrorCode
  objective(Tao, Vec x, PetscReal *value, void *)
  {
    // The objective function is 1^T x = sum_i x_i
    PetscScalar sum;
    PetscFunctionBegin;
    PetscCall(VecSum(x, &sum));
    *value = PetscRealPart(sum);
    PetscFunctionReturn(PETSC_SUCCESS);
  }



  PetscErrorCode
  gradient(Tao, Vec, Vec gradient_vector, void *)
  {
    // The objective function is 1^T x, so the gradient is the vector of ones
    PetscFunctionBegin;
    PetscCall(VecSet(gradient_vector, 1.));
    PetscFunctionReturn(PETSC_SUCCESS);
  }



  PetscErrorCode
  hessian(Tao, Vec, Mat hessian_matrix, Mat, void *)
  {
    // It's a linear program, so the Hessian is the zero matrix
    PetscFunctionBegin;
    PetscCall(MatZeroEntries(hessian_matrix));
    PetscFunctionReturn(PETSC_SUCCESS);
  }



  PetscErrorCode
  equality_constraints(Tao, Vec, Vec, void *)
  {
    // We have no equality constraints
    return PETSC_SUCCESS;
  }



  PetscErrorCode
  equality_jacobian(Tao, Vec, Mat jacobian, Mat, void *)
  {
    // We have no equality constraints, so the Jacobian is an empty matrix
    PetscFunctionBegin;
    PetscCall(MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY));
    PetscFunctionReturn(PETSC_SUCCESS);
  }



  PetscErrorCode
  inequality_constraints(Tao, Vec x, Vec values, void *)
  {
    // We need to compute h(x)_j = x_{j+1}-x_j-1, which needs to be >=0
    const PetscScalar *x_entries;
    PetscScalar       *values_entries;
    PetscFunctionBegin;
    PetscCall(VecGetArrayRead(x, &x_entries));
    PetscCall(VecGetArray(values, &values_entries));
    for (PetscInt index = 0; index < n_inequalities; ++index)
      values_entries[index] = x_entries[index + 1] - x_entries[index] - 1.;
    PetscCall(VecRestoreArray(values, &values_entries));
    PetscCall(VecRestoreArrayRead(x, &x_entries));
    PetscFunctionReturn(PETSC_SUCCESS);
  }



  PetscErrorCode
  inequality_jacobian(Tao, Vec, Mat jacobian, Mat, void *)
  {
    // We have h(x)_j = x_{j+1}-x_j-1, so the Jacobian J_{jk} is
    // J_{j,j+1}=1, J_{jj}=-1, and zero everywhere else
    PetscFunctionBegin;
    PetscCall(MatZeroEntries(jacobian));
    for (PetscInt index = 0; index < n_inequalities; ++index)
      {
        PetscCall(MatSetValue(jacobian, index, index, -1., INSERT_VALUES));
        PetscCall(MatSetValue(jacobian, index, index + 1, 1., INSERT_VALUES));
      }
    PetscCall(MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
} // namespace



int
main(int argc, char **argv)
{
  const PetscErrorCode initialization_error =
    PetscInitialize(&argc, &argv, nullptr, nullptr);
  if (initialization_error != PETSC_SUCCESS)
    return initialization_error;
  PetscCallAbort(PETSC_COMM_WORLD, PetscSetFPTrap(PETSC_FP_TRAP_ON));

  // Set up matrix and vector objects:
  Vec x;
  Vec gradient_vector;
  Vec equality_constraint_values;
  Vec inequality_constraint_values;
  Vec lower_bounds;
  Vec upper_bounds;
  Mat hessian_matrix;
  Mat equality_jacobian_matrix;
  Mat inequality_jacobian_matrix;

  Tao tao;
  PetscCallAbort(MPI_COMM_SELF, VecCreateSeq(MPI_COMM_SELF, n_variables, &x));
  PetscCallAbort(MPI_COMM_SELF,
                 VecCreateSeq(MPI_COMM_SELF, n_variables, &gradient_vector));
  PetscCallAbort(MPI_COMM_SELF,
                 VecCreateSeq(MPI_COMM_SELF, 0, &equality_constraint_values));
  PetscCallAbort(MPI_COMM_SELF,
                 VecCreateSeq(MPI_COMM_SELF,
                              n_inequalities,
                              &inequality_constraint_values));
  PetscCallAbort(MPI_COMM_SELF,
                 VecCreateSeq(MPI_COMM_SELF, n_variables, &lower_bounds));
  PetscCallAbort(MPI_COMM_SELF,
                 VecCreateSeq(MPI_COMM_SELF, n_variables, &upper_bounds));

  // Initialize the solution vector with x_i=i+10+U(0,100), which is a feasible
  // starting point. Lower bounds are x_i>=0 and we have upper bounds
  // x_i<=infinity.
  for (PetscInt index = 0; index < n_variables; ++index)
    {
      std::mt19937                           generator(index);
      std::uniform_real_distribution<double> distribution(0., 100.);
      PetscCallAbort(MPI_COMM_SELF,
                     VecSetValue(x,
                                 index,
                                 index + 10. + distribution(generator),
                                 INSERT_VALUES));
      PetscCallAbort(MPI_COMM_SELF,
                     VecSetValue(lower_bounds, index, 0., INSERT_VALUES));
      PetscCallAbort(
        MPI_COMM_SELF,
        VecSetValue(upper_bounds, index, PETSC_INFINITY, INSERT_VALUES));
    }
  PetscCallAbort(MPI_COMM_SELF, VecAssemblyBegin(x));
  PetscCallAbort(MPI_COMM_SELF, VecAssemblyEnd(x));
  PetscCallAbort(MPI_COMM_SELF, VecAssemblyBegin(lower_bounds));
  PetscCallAbort(MPI_COMM_SELF, VecAssemblyEnd(lower_bounds));
  PetscCallAbort(MPI_COMM_SELF, VecAssemblyBegin(upper_bounds));
  PetscCallAbort(MPI_COMM_SELF, VecAssemblyEnd(upper_bounds));

  // Create the Hessian as a diagonal zero matrix
  PetscCallAbort(
    MPI_COMM_SELF,
    MatCreateSeqAIJ(
      MPI_COMM_SELF, n_variables, n_variables, 1, nullptr, &hessian_matrix));
  for (PetscInt index = 0; index < n_variables; ++index)
    PetscCallAbort(
      MPI_COMM_SELF,
      MatSetValue(hessian_matrix, index, index, 0., INSERT_VALUES));
  PetscCallAbort(MPI_COMM_SELF,
                 MatAssemblyBegin(hessian_matrix, MAT_FINAL_ASSEMBLY));
  PetscCallAbort(MPI_COMM_SELF,
                 MatAssemblyEnd(hessian_matrix, MAT_FINAL_ASSEMBLY));

  // Create an empty equality Jacobian matrix (we have no equality constraints)
  PetscCallAbort(
    MPI_COMM_SELF,
    MatCreateSeqAIJ(
      MPI_COMM_SELF, 0, n_variables, 0, nullptr, &equality_jacobian_matrix));
  PetscCallAbort(MPI_COMM_SELF,
                 MatAssemblyBegin(equality_jacobian_matrix,
                                  MAT_FINAL_ASSEMBLY));
  PetscCallAbort(MPI_COMM_SELF,
                 MatAssemblyEnd(equality_jacobian_matrix, MAT_FINAL_ASSEMBLY));

  // For inequalities, we have one inequality for each of x_2...x_n,
  // with 2 nonzeros per row in the matrix:
  PetscCallAbort(MPI_COMM_SELF,
                 MatCreateSeqAIJ(MPI_COMM_SELF,
                                 n_inequalities,
                                 n_variables,
                                 2,
                                 nullptr,
                                 &inequality_jacobian_matrix));

  // Now set up the TAO object and tie it to the various objects above:
  PetscCallAbort(MPI_COMM_SELF, TaoCreate(MPI_COMM_SELF, &tao));
  PetscCallAbort(MPI_COMM_SELF, TaoSetType(tao, TAOALMM));
  PetscCallAbort(MPI_COMM_SELF, TaoALMMSetType(tao, TAO_ALMM_PHR));
  PetscCallAbort(MPI_COMM_SELF,
                 PetscOptionsSetValue(nullptr,
                                      "-tao_almm_subsolver_tao_type",
                                      "bncg"));
  PetscCallAbort(MPI_COMM_SELF,
                 PetscOptionsSetValue(nullptr, "-tao_almm_mu_init", "1."));
  PetscCallAbort(MPI_COMM_SELF,
                 PetscOptionsSetValue(nullptr, "-tao_almm_mu_factor", "2."));
#if PETSC_VERSION_GE(3, 17, 0)
  PetscCallAbort(MPI_COMM_SELF, TaoSetSolution(tao, x));
#else
  PetscCallAbort(MPI_COMM_SELF, TaoSetInitialVector(tao, x));
#endif
  PetscCallAbort(MPI_COMM_SELF,
                 TaoSetVariableBounds(tao, lower_bounds, upper_bounds));
#if PETSC_VERSION_GE(3, 17, 0)
  PetscCallAbort(MPI_COMM_SELF, TaoSetObjective(tao, objective, nullptr));
  PetscCallAbort(MPI_COMM_SELF,
                 TaoSetGradient(tao, gradient_vector, gradient, nullptr));
#else
  PetscCallAbort(MPI_COMM_SELF,
                 TaoSetObjectiveRoutine(tao, objective, nullptr));
  PetscCallAbort(MPI_COMM_SELF, TaoSetGradientRoutine(tao, gradient, nullptr));
#endif
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
  PetscCallAbort(
    MPI_COMM_SELF,
    TaoSetInequalityConstraintsRoutine(
      tao, inequality_constraint_values, inequality_constraints, nullptr));
  PetscCallAbort(MPI_COMM_SELF,
                 TaoSetJacobianInequalityRoutine(tao,
                                                 inequality_jacobian_matrix,
                                                 inequality_jacobian_matrix,
                                                 inequality_jacobian,
                                                 nullptr));
#if PETSC_VERSION_GE(3, 17, 0)
  PetscCallAbort(
    MPI_COMM_SELF,
    TaoSetHessian(tao, hessian_matrix, hessian_matrix, hessian, nullptr));
#else
  PetscCallAbort(MPI_COMM_SELF,
                 TaoSetHessianRoutine(
                   tao, hessian_matrix, hessian_matrix, hessian, nullptr));
#endif
  PetscCallAbort(MPI_COMM_SELF, TaoSetFromOptions(tao));

  // Then solve the problem and output how many iterations that took:
  PetscCallAbort(MPI_COMM_SELF, TaoSolve(tao));

  PetscInt iterations;
  PetscCallAbort(MPI_COMM_SELF, TaoGetIterationNumber(tao, &iterations));
  std::cout << "The optimizer found the solution in " << iterations
            << " iterations." << std::endl;

  // We know what the solution should be (x_i=i). Print it, then assert that it
  // is correct:
  const PetscScalar *solution;
  PetscCallAbort(MPI_COMM_SELF, VecGetArrayRead(x, &solution));
  for (PetscInt index = 0; index < n_variables; ++index)
    std::cout << "x[" << index << "]=" << solution[index]
              << ", error=" << (index - solution[index]) << std::endl;
  for (PetscInt index = 0; index < n_variables; ++index)
    PetscCheckAbort(std::abs(PetscRealPart(solution[index]) - index) < 1e-6,
                    MPI_COMM_SELF,
                    PETSC_ERR_PLIB,
                    "Solution component %" PetscInt_FMT
                    " does not equal %" PetscInt_FMT ".",
                    index,
                    index);
  PetscCallAbort(MPI_COMM_SELF, VecRestoreArrayRead(x, &solution));

  // Final clean-ups:
  PetscCallAbort(MPI_COMM_SELF, TaoDestroy(&tao));
  PetscCallAbort(MPI_COMM_SELF, MatDestroy(&inequality_jacobian_matrix));
  PetscCallAbort(MPI_COMM_SELF, MatDestroy(&equality_jacobian_matrix));
  PetscCallAbort(MPI_COMM_SELF, MatDestroy(&hessian_matrix));
  PetscCallAbort(MPI_COMM_SELF, VecDestroy(&upper_bounds));
  PetscCallAbort(MPI_COMM_SELF, VecDestroy(&lower_bounds));
  PetscCallAbort(MPI_COMM_SELF, VecDestroy(&inequality_constraint_values));
  PetscCallAbort(MPI_COMM_SELF, VecDestroy(&equality_constraint_values));
  PetscCallAbort(MPI_COMM_SELF, VecDestroy(&gradient_vector));
  PetscCallAbort(MPI_COMM_SELF, VecDestroy(&x));
  PetscCallAbort(PETSC_COMM_WORLD, PetscFinalize());
}
