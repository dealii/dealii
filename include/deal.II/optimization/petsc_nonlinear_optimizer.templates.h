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

#ifndef dealii_petsc_nonlinear_optimizer_templates_h
#define dealii_petsc_nonlinear_optimizer_templates_h

#include <deal.II/base/config.h>

#if defined(DEAL_II_WITH_PETSC)
#  if DEAL_II_PETSC_VERSION_GTE(3, 13, 0)
#    include <deal.II/base/scope_exit.h>

#    include <deal.II/lac/exceptions.h>

#    include <deal.II/optimization/petsc_nonlinear_optimizer.h>

#    include <utility>
#  endif
#endif

DEAL_II_NAMESPACE_OPEN

#if defined(DEAL_II_WITH_PETSC)
#  if DEAL_II_PETSC_VERSION_GTE(3, 13, 0)
namespace PETScWrappers
{
  namespace internal
  {
    /**
     * Assert that a PETSc call made by the nonlinear optimizer succeeded.
     */
    inline void
    assert_nonlinear_optimizer_petsc(const PetscErrorCode error_code)
    {
      AssertThrow(error_code == 0, ExcPETScError(error_code));
    }



    /**
     * Return success for an empty nonlinear optimizer constraint set.
     */
    inline PetscErrorCode
    empty_nonlinear_optimizer_constraints(Tao, Vec, Vec, void *)
    {
      return 0;
    }



    /**
     * Assemble the Jacobian associated with an empty constraint set.
     */
    inline PetscErrorCode
    empty_nonlinear_optimizer_jacobian(Tao, Vec, Mat jacobian, Mat, void *)
    {
      PetscErrorCode error_code =
        MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY);
      if (error_code != 0)
        return error_code;
      error_code = MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY);
      return error_code;
    }



    /**
     * Invoke a nonlinear optimizer callback and defer any exception to the
     * caller through the pending-exception argument.
     */
    template <typename Function>
    PetscErrorCode
    call_and_capture_nonlinear_optimizer_exception(
      std::exception_ptr &pending_exception,
      const Function     &function)
    {
      try
        {
          function();
          return 0;
        }
      catch (...)
        {
          pending_exception = std::current_exception();
          return PETSC_ERR_USER;
        }
    }
  } // namespace internal



  template <typename VectorType, typename MatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
  NonlinearOptimizer<VectorType, MatrixType>::NonlinearOptimizer(
    const AdditionalData &data,
    const MPI_Comm        mpi_comm)
    : tao(nullptr)
  {
    internal::assert_nonlinear_optimizer_petsc(TaoCreate(mpi_comm, &tao));
    internal::assert_nonlinear_optimizer_petsc(
      TaoSetApplicationContext(tao, this));
    internal::assert_nonlinear_optimizer_petsc(TaoSetType(tao, TAOALMM));
    Tao almm_subsolver;
    internal::assert_nonlinear_optimizer_petsc(
      TaoALMMGetSubsolver(tao, &almm_subsolver));
    switch (data.almm_subsolver)
      {
        case AdditionalData::ALMMSubsolver::bqnls:
          internal::assert_nonlinear_optimizer_petsc(
            TaoSetType(almm_subsolver, TAOBQNLS));
          break;

        case AdditionalData::ALMMSubsolver::bncg:
          internal::assert_nonlinear_optimizer_petsc(
            TaoSetType(almm_subsolver, TAOBNCG));
          break;
      }

    if (!data.options_prefix.empty())
      internal::assert_nonlinear_optimizer_petsc(
        TaoSetOptionsPrefix(tao, data.options_prefix.c_str()));

    const PetscReal gradient_tolerance =
      data.gradient_tolerance > 0. ? data.gradient_tolerance : PETSC_DEFAULT;
    internal::assert_nonlinear_optimizer_petsc(
      TaoSetTolerances(tao, gradient_tolerance, PETSC_DEFAULT, PETSC_DEFAULT));
    if (data.constraint_tolerance > 0.)
      internal::assert_nonlinear_optimizer_petsc(TaoSetConstraintTolerances(
        tao, data.constraint_tolerance, PETSC_DEFAULT));
    if (data.maximum_iterations >= 0)
      internal::assert_nonlinear_optimizer_petsc(
        TaoSetMaximumIterations(tao, data.maximum_iterations));
  }



  template <typename VectorType, typename MatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
  NonlinearOptimizer<VectorType, MatrixType>::~NonlinearOptimizer()
  {
    internal::assert_nonlinear_optimizer_petsc(TaoDestroy(&tao));
  }



  template <typename VectorType, typename MatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
  unsigned int NonlinearOptimizer<VectorType, MatrixType>::solve(
    VectorType                       &solution,
    MatrixType                       &hessian_matrix,
    const std::optional<Constraints> &equality,
    const std::optional<Constraints> &inequality,
    const std::optional<Bounds>      &bounds)
  {
    AssertThrow(objective, ExcMessage("An objective callback is required."));
    AssertThrow(gradient, ExcMessage("A gradient callback is required."));
    AssertThrow(hessian, ExcMessage("A Hessian callback is required."));

    // Error checking: If bounds are supplied, check that the initial
    // guess satisfies them:
    if (bounds.has_value())
      {
        Vec bound_distance = nullptr;
        internal::assert_nonlinear_optimizer_petsc(
          VecDuplicate(solution.petsc_vector(), &bound_distance));
        ScopeExit destroy_bound_distance([&]() {
          internal::assert_nonlinear_optimizer_petsc(
            VecDestroy(&bound_distance));
        });

        PetscReal minimum_lower_bound_distance;
        internal::assert_nonlinear_optimizer_petsc(
          VecWAXPY(bound_distance,
                   -1.,
                   static_cast<const Vec &>(bounds->lower),
                   solution.petsc_vector()));
        internal::assert_nonlinear_optimizer_petsc(
          VecMin(bound_distance, nullptr, &minimum_lower_bound_distance));
        AssertThrow(minimum_lower_bound_distance >= 0.,
                    ExcMessage(
                      "The initial guess violates a lower variable bound."));

        PetscReal minimum_upper_bound_distance;
        internal::assert_nonlinear_optimizer_petsc(
          VecWAXPY(bound_distance,
                   -1.,
                   solution.petsc_vector(),
                   static_cast<const Vec &>(bounds->upper)));
        internal::assert_nonlinear_optimizer_petsc(
          VecMin(bound_distance, nullptr, &minimum_upper_bound_distance));
        AssertThrow(minimum_upper_bound_distance >= 0.,
                    ExcMessage(
                      "The initial guess violates an upper variable bound."));
      }

#    if DEAL_II_PETSC_VERSION_GTE(3, 17, 0)
    internal::assert_nonlinear_optimizer_petsc(
      TaoSetSolution(tao, solution.petsc_vector()));
#    else
    internal::assert_nonlinear_optimizer_petsc(
      TaoSetInitialVector(tao, solution.petsc_vector()));
#    endif

    Vec       dummy_equality_values     = nullptr;
    Mat       dummy_equality_jacobian   = nullptr;
    Vec       dummy_inequality_values   = nullptr;
    Mat       dummy_inequality_jacobian = nullptr;
    Vec       dummy_lower_bounds        = nullptr;
    Vec       dummy_upper_bounds        = nullptr;
    ScopeExit destroy_dummy_objects([&]() {
      internal::assert_nonlinear_optimizer_petsc(
        VecDestroy(&dummy_equality_values));
      internal::assert_nonlinear_optimizer_petsc(
        MatDestroy(&dummy_equality_jacobian));
      internal::assert_nonlinear_optimizer_petsc(
        VecDestroy(&dummy_inequality_values));
      internal::assert_nonlinear_optimizer_petsc(
        MatDestroy(&dummy_inequality_jacobian));
      internal::assert_nonlinear_optimizer_petsc(
        VecDestroy(&dummy_lower_bounds));
      internal::assert_nonlinear_optimizer_petsc(
        VecDestroy(&dummy_upper_bounds));
    });

    // If equality or inequality constraints are supplied, let TAO know that
    // they are present and provide the appropriate callbacks.
    if (equality.has_value())
      {
        AssertThrow(equality_constraints && equality_jacobian,
                    ExcMessage(
                      "If you provide an equality constraint object, you must "
                      "also set the equality constraint callbacks."));
        internal::assert_nonlinear_optimizer_petsc(
          TaoSetEqualityConstraintsRoutine(
            tao,
            const_cast<VectorType &>(equality->values).petsc_vector(),
            &NonlinearOptimizer::equality_constraints_callback,
            this));
        internal::assert_nonlinear_optimizer_petsc(
          TaoSetJacobianEqualityRoutine(
            tao,
            const_cast<MatrixType &>(equality->jacobian).petsc_matrix(),
            const_cast<MatrixType &>(equality->jacobian).petsc_matrix(),
            &NonlinearOptimizer::equality_jacobian_callback,
            this));
      }

    if (inequality.has_value())
      {
        AssertThrow(
          inequality_constraints && inequality_jacobian,
          ExcMessage("If you provide an inequality constraint object, you must "
                     "also set the inequality constraint callbacks."));
        internal::assert_nonlinear_optimizer_petsc(
          TaoSetInequalityConstraintsRoutine(
            tao,
            const_cast<VectorType &>(inequality->values).petsc_vector(),
            &NonlinearOptimizer::inequality_constraints_callback,
            this));
        internal::assert_nonlinear_optimizer_petsc(
          TaoSetJacobianInequalityRoutine(
            tao,
            const_cast<MatrixType &>(inequality->jacobian).petsc_matrix(),
            const_cast<MatrixType &>(inequality->jacobian).petsc_matrix(),
            &NonlinearOptimizer::inequality_jacobian_callback,
            this));
      }

    // TAOALMM requires constraint objects, even for an unconstrained problem.
    // When neither constraint type is supplied, use empty equality and
    // inequality constraints to represent that case.
    if (!equality.has_value() && !inequality.has_value())
      {
        const Vec solution_vector = solution.petsc_vector();
        PetscInt  local_size;
        PetscInt  global_size;
        internal::assert_nonlinear_optimizer_petsc(
          VecGetLocalSize(solution_vector, &local_size));
        internal::assert_nonlinear_optimizer_petsc(
          VecGetSize(solution_vector, &global_size));

        const MPI_Comm mpi_comm = solution.get_mpi_communicator();
        internal::assert_nonlinear_optimizer_petsc(
          VecCreateMPI(mpi_comm, 0, 0, &dummy_equality_values));
        internal::assert_nonlinear_optimizer_petsc(
          MatCreateAIJ(mpi_comm,
                       0,
                       local_size,
                       0,
                       global_size,
                       0,
                       nullptr,
                       0,
                       nullptr,
                       &dummy_equality_jacobian));
        internal::assert_nonlinear_optimizer_petsc(
          MatAssemblyBegin(dummy_equality_jacobian, MAT_FINAL_ASSEMBLY));
        internal::assert_nonlinear_optimizer_petsc(
          MatAssemblyEnd(dummy_equality_jacobian, MAT_FINAL_ASSEMBLY));
        internal::assert_nonlinear_optimizer_petsc(
          TaoSetEqualityConstraintsRoutine(
            tao,
            dummy_equality_values,
            &internal::empty_nonlinear_optimizer_constraints,
            nullptr));
        internal::assert_nonlinear_optimizer_petsc(
          TaoSetJacobianEqualityRoutine(
            tao,
            dummy_equality_jacobian,
            dummy_equality_jacobian,
            &internal::empty_nonlinear_optimizer_jacobian,
            nullptr));

        internal::assert_nonlinear_optimizer_petsc(
          VecCreateMPI(mpi_comm, 0, 0, &dummy_inequality_values));
        internal::assert_nonlinear_optimizer_petsc(
          MatCreateAIJ(mpi_comm,
                       0,
                       local_size,
                       0,
                       global_size,
                       0,
                       nullptr,
                       0,
                       nullptr,
                       &dummy_inequality_jacobian));
        internal::assert_nonlinear_optimizer_petsc(
          MatAssemblyBegin(dummy_inequality_jacobian, MAT_FINAL_ASSEMBLY));
        internal::assert_nonlinear_optimizer_petsc(
          MatAssemblyEnd(dummy_inequality_jacobian, MAT_FINAL_ASSEMBLY));
        internal::assert_nonlinear_optimizer_petsc(
          TaoSetInequalityConstraintsRoutine(
            tao,
            dummy_inequality_values,
            &internal::empty_nonlinear_optimizer_constraints,
            nullptr));
        internal::assert_nonlinear_optimizer_petsc(
          TaoSetJacobianInequalityRoutine(
            tao,
            dummy_inequality_jacobian,
            dummy_inequality_jacobian,
            &internal::empty_nonlinear_optimizer_jacobian,
            nullptr));
      }

    // If present, set the variable bounds. If not present, create dummy PETSc
    // objects to satisfy TAO's requirement that bounds be present. These
    // objects will not be used, since the callbacks will be set to empty
    // functions.
    if (bounds.has_value())
      internal::assert_nonlinear_optimizer_petsc(
        TaoSetVariableBounds(tao,
                             static_cast<const Vec &>(bounds->lower),
                             static_cast<const Vec &>(bounds->upper)));
    else
      {
        const Vec solution_vector = solution.petsc_vector();
        internal::assert_nonlinear_optimizer_petsc(
          VecDuplicate(solution_vector, &dummy_lower_bounds));
        internal::assert_nonlinear_optimizer_petsc(
          VecDuplicate(solution_vector, &dummy_upper_bounds));
        internal::assert_nonlinear_optimizer_petsc(
          VecSet(dummy_lower_bounds, PETSC_NINFINITY));
        internal::assert_nonlinear_optimizer_petsc(
          VecSet(dummy_upper_bounds, PETSC_INFINITY));
        internal::assert_nonlinear_optimizer_petsc(
          TaoSetVariableBounds(tao, dummy_lower_bounds, dummy_upper_bounds));
      }

      // Other setup:
#    if DEAL_II_PETSC_VERSION_GTE(3, 17, 0)
    internal::assert_nonlinear_optimizer_petsc(
      TaoSetObjective(tao, &NonlinearOptimizer::objective_callback, this));
    internal::assert_nonlinear_optimizer_petsc(TaoSetGradient(
      tao, nullptr, &NonlinearOptimizer::gradient_callback, this));
    internal::assert_nonlinear_optimizer_petsc(
      TaoSetHessian(tao,
                    hessian_matrix.petsc_matrix(),
                    hessian_matrix.petsc_matrix(),
                    &NonlinearOptimizer::hessian_callback,
                    this));
#    else
    internal::assert_nonlinear_optimizer_petsc(TaoSetObjectiveRoutine(
      tao, &NonlinearOptimizer::objective_callback, this));
    internal::assert_nonlinear_optimizer_petsc(
      TaoSetGradientRoutine(tao, &NonlinearOptimizer::gradient_callback, this));
    internal::assert_nonlinear_optimizer_petsc(
      TaoSetHessianRoutine(tao,
                           hessian_matrix.petsc_matrix(),
                           hessian_matrix.petsc_matrix(),
                           &NonlinearOptimizer::hessian_callback,
                           this));
#    endif
    internal::assert_nonlinear_optimizer_petsc(TaoSetFromOptions(tao));

    // Solve the optimization problem. If the user callback throws an exception,
    // it will be captured and rethrown here.
    internal::assert_nonlinear_optimizer_petsc(TaoSolve(tao));

    if (pending_exception != nullptr)
      std::rethrow_exception(std::exchange(pending_exception, nullptr));

    // Obtain the number of iterations taken by the solver and return it to the
    // caller:
    PetscInt iterations;
    internal::assert_nonlinear_optimizer_petsc(
      TaoGetIterationNumber(tao, &iterations));
    return iterations;
  }



  template <typename VectorType, typename MatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
  Tao NonlinearOptimizer<VectorType, MatrixType>::petsc_tao()
  {
    return tao;
  }



  template <typename VectorType, typename MatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
  PetscErrorCode NonlinearOptimizer<VectorType, MatrixType>::objective_callback(
    Tao,
    Vec        x,
    PetscReal *value,
    void      *context)
  {
    auto &optimizer = *static_cast<NonlinearOptimizer *>(context);
    return internal::call_and_capture_nonlinear_optimizer_exception(
      optimizer.pending_exception, [&]() {
        VectorType x_wrapper(x);
        *value = optimizer.objective(x_wrapper);
      });
  }



  template <typename VectorType, typename MatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
  PetscErrorCode
    NonlinearOptimizer<VectorType, MatrixType>::gradient_callback(Tao,
                                                                  Vec x,
                                                                  Vec gradient,
                                                                  void *context)
  {
    auto &optimizer = *static_cast<NonlinearOptimizer *>(context);
    return internal::call_and_capture_nonlinear_optimizer_exception(
      optimizer.pending_exception, [&]() {
        VectorType x_wrapper(x);
        VectorType gradient_wrapper(gradient);
        optimizer.gradient(x_wrapper, gradient_wrapper);
      });
  }



  template <typename VectorType, typename MatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
  PetscErrorCode NonlinearOptimizer<VectorType, MatrixType>::hessian_callback(
    Tao,
    Vec   x,
    Mat   hessian,
    Mat   preconditioner,
    void *context)
  {
    auto &optimizer = *static_cast<NonlinearOptimizer *>(context);
    return internal::call_and_capture_nonlinear_optimizer_exception(
      optimizer.pending_exception, [&]() {
        VectorType x_wrapper(x);
        MatrixType hessian_wrapper(hessian);
        MatrixType preconditioner_wrapper(preconditioner);
        optimizer.hessian(x_wrapper, hessian_wrapper, preconditioner_wrapper);
      });
  }



  template <typename VectorType, typename MatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
  PetscErrorCode
    NonlinearOptimizer<VectorType, MatrixType>::equality_constraints_callback(
      Tao,
      Vec   x,
      Vec   values,
      void *context)
  {
    auto &optimizer = *static_cast<NonlinearOptimizer *>(context);
    return internal::call_and_capture_nonlinear_optimizer_exception(
      optimizer.pending_exception, [&]() {
        VectorType x_wrapper(x);
        VectorType values_wrapper(values);
        optimizer.equality_constraints(x_wrapper, values_wrapper);
      });
  }



  template <typename VectorType, typename MatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
  PetscErrorCode
    NonlinearOptimizer<VectorType, MatrixType>::equality_jacobian_callback(
      Tao,
      Vec   x,
      Mat   jacobian,
      Mat   preconditioner,
      void *context)
  {
    auto &optimizer = *static_cast<NonlinearOptimizer *>(context);
    return internal::call_and_capture_nonlinear_optimizer_exception(
      optimizer.pending_exception, [&]() {
        VectorType x_wrapper(x);
        MatrixType jacobian_wrapper(jacobian);
        MatrixType preconditioner_wrapper(preconditioner);
        optimizer.equality_jacobian(x_wrapper,
                                    jacobian_wrapper,
                                    preconditioner_wrapper);
      });
  }



  template <typename VectorType, typename MatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
  PetscErrorCode
    NonlinearOptimizer<VectorType, MatrixType>::inequality_constraints_callback(
      Tao,
      Vec   x,
      Vec   values,
      void *context)
  {
    auto &optimizer = *static_cast<NonlinearOptimizer *>(context);
    return internal::call_and_capture_nonlinear_optimizer_exception(
      optimizer.pending_exception, [&]() {
        VectorType x_wrapper(x);
        VectorType values_wrapper(values);
        optimizer.inequality_constraints(x_wrapper, values_wrapper);
      });
  }



  template <typename VectorType, typename MatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
  PetscErrorCode
    NonlinearOptimizer<VectorType, MatrixType>::inequality_jacobian_callback(
      Tao,
      Vec   x,
      Mat   jacobian,
      Mat   preconditioner,
      void *context)
  {
    auto &optimizer = *static_cast<NonlinearOptimizer *>(context);
    return internal::call_and_capture_nonlinear_optimizer_exception(
      optimizer.pending_exception, [&]() {
        VectorType x_wrapper(x);
        MatrixType jacobian_wrapper(jacobian);
        MatrixType preconditioner_wrapper(preconditioner);
        optimizer.inequality_jacobian(x_wrapper,
                                      jacobian_wrapper,
                                      preconditioner_wrapper);
      });
  }
} // namespace PETScWrappers
#  endif
#endif

DEAL_II_NAMESPACE_CLOSE

#endif
