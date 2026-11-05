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

#ifndef dealii_petsc_linear_optimizer_templates_h
#define dealii_petsc_linear_optimizer_templates_h

#include <deal.II/base/config.h>

#if defined(DEAL_II_WITH_PETSC)
#  if DEAL_II_PETSC_VERSION_GTE(3, 13, 0)
#    include <deal.II/base/scope_exit.h>

#    include <deal.II/optimization/petsc_linear_optimizer.h>
#    include <deal.II/optimization/petsc_nonlinear_optimizer.templates.h>

#    include <optional>
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
     * Prefer a MatrixType-specific reinit(Mat) member function when one is
     * available. The trailing return type is void, while its decltype
     * expression uses SFINAE to remove this overload if that member function
     * does not exist; passing 0 then permits the long overload below to be
     * selected as a fallback.
     */
    template <typename MatrixType>
    auto
    reinit_linear_program_matrix(MatrixType &matrix, Mat mat, int)
      -> decltype(matrix.reinit(mat), void())
    {
      matrix.reinit(mat);
    }



    /**
     * Fall back to the MatrixBase interface when MatrixType does not expose
     * its own reinit(Mat) member function.
     */
    template <typename MatrixType>
    void
    reinit_linear_program_matrix(MatrixType &matrix, Mat mat, long)
    {
      static_cast<MatrixBase &>(matrix).reinit(mat);
    }
  } // namespace internal



  template <typename VectorType, typename MatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
  LinearOptimizer<VectorType, MatrixType>::LinearOptimizer(
    const AdditionalData &data,
    const MPI_Comm        mpi_comm)
    : optimizer(data, mpi_comm)
  {}



  template <typename VectorType, typename MatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
  unsigned int LinearOptimizer<VectorType, MatrixType>::solve(
    VectorType                              &solution,
    const VectorType                        &objective,
    const ObserverPointer<MatrixType>       &equality_matrix,
    const ObserverPointer<const VectorType> &equality_right_hand_side,
    const ObserverPointer<MatrixType>       &inequality_matrix,
    const ObserverPointer<const VectorType> &inequality_right_hand_side,
    const std::optional<
      typename NonlinearOptimizer<VectorType, MatrixType>::Bounds> &bounds)
  {
    AssertThrow(
      (equality_matrix == nullptr) == (equality_right_hand_side == nullptr),
      ExcMessage(
        "Equality matrix and right-hand side must both either be set or absent."));
    AssertThrow(
      (inequality_matrix == nullptr) == (inequality_right_hand_side == nullptr),
      ExcMessage(
        "Inequality matrix and right-hand side must both either be set or absent."));

    // Create a Hessian matrix to pass to the nonlinear optimizer. This is
    // required because the nonlinear optimizer interface requires a Hessian
    // matrix, even though it is not used for linear problems.
    const PetscInt local_size  = solution.locally_owned_size();
    const PetscInt global_size = solution.size();
    Mat            hessian     = nullptr;
    ScopeExit      cleanup_hessian([&]() {
      internal::assert_nonlinear_optimizer_petsc(MatDestroy(&hessian));
    });
    internal::assert_nonlinear_optimizer_petsc(
      MatCreateAIJ(solution.get_mpi_communicator(),
                   local_size,
                   local_size,
                   global_size,
                   global_size,
                   1,
                   nullptr,
                   1,
                   nullptr,
                   &hessian));
    internal::assert_nonlinear_optimizer_petsc(MatSetUp(hessian));

    // Zero out the Hessian matrix, since it is not used for linear problems.
    // One could imagine that a MatShell could be used here, implementing
    // a zero matrix, but that run afoul of TAO's requirement that the Hessian
    // be an "assembled" matrix, which is not the case for a MatShell.
    PetscInt first_owned_row;
    PetscInt last_owned_row;
    internal::assert_nonlinear_optimizer_petsc(
      MatGetOwnershipRange(hessian, &first_owned_row, &last_owned_row));
    const PetscScalar zero = 0.;
    for (PetscInt row = first_owned_row; row < last_owned_row; ++row)
      internal::assert_nonlinear_optimizer_petsc(
        MatSetValues(hessian, 1, &row, 1, &row, &zero, INSERT_VALUES));
    internal::assert_nonlinear_optimizer_petsc(
      MatAssemblyBegin(hessian, MAT_FINAL_ASSEMBLY));
    internal::assert_nonlinear_optimizer_petsc(
      MatAssemblyEnd(hessian, MAT_FINAL_ASSEMBLY));

    MatrixType hessian_wrapper(hessian);

    // Set up the nonlinear optimizer with callbacks that implement the linear
    // problem. The objective and gradient are linear, and the Hessian is zero.
    const ObserverPointer<const VectorType> objective_pointer(&objective);
    optimizer.objective =
      [objective_pointer](const VectorType &x) -> PetscReal {
      PetscScalar objective_value;
      internal::assert_nonlinear_optimizer_petsc(
        VecDot(static_cast<const Vec &>(x),
               static_cast<const Vec &>(*objective_pointer),
               &objective_value));
      return PetscRealPart(objective_value);
    };
    optimizer.gradient = [objective_pointer](const VectorType &,
                                             VectorType &gradient) {
      internal::assert_nonlinear_optimizer_petsc(
        VecCopy(static_cast<const Vec &>(*objective_pointer),
                gradient.petsc_vector()));
    };
    optimizer.hessian =
      [](const VectorType &, MatrixType &hessian, MatrixType &) {
        internal::assert_nonlinear_optimizer_petsc(
          MatZeroEntries(hessian.petsc_matrix()));
      };

    // Translate the equality and inequality constraints into the nonlinear
    // optimizer's interface. This requires creating copies of the constraint
    // matrices.
    std::optional<
      typename NonlinearOptimizer<VectorType, MatrixType>::Constraints>
      equality;
    if (equality_matrix != nullptr)
      {
        equality.emplace();
        equality->values = *equality_right_hand_side;
        Mat equality_jacobian;
        internal::assert_nonlinear_optimizer_petsc(
          MatDuplicate(equality_matrix->petsc_matrix(),
                       MAT_COPY_VALUES,
                       &equality_jacobian));
        internal::reinit_linear_program_matrix(equality->jacobian,
                                               equality_jacobian,
                                               0);
        internal::assert_nonlinear_optimizer_petsc(
          MatDestroy(&equality_jacobian));
        optimizer.equality_constraints =
          [equality_matrix, equality_right_hand_side](const VectorType &x,
                                                      VectorType &values) {
            internal::assert_nonlinear_optimizer_petsc(
              MatMult(equality_matrix->petsc_matrix(),
                      static_cast<const Vec &>(x),
                      values.petsc_vector()));
            internal::assert_nonlinear_optimizer_petsc(
              VecAXPY(values.petsc_vector(),
                      -1.,
                      static_cast<const Vec &>(*equality_right_hand_side)));
          };
        optimizer.equality_jacobian = [equality_matrix](const VectorType &,
                                                        MatrixType &jacobian,
                                                        MatrixType &) {
          internal::assert_nonlinear_optimizer_petsc(
            MatCopy(equality_matrix->petsc_matrix(),
                    jacobian.petsc_matrix(),
                    SAME_NONZERO_PATTERN));
        };
      }

    std::optional<
      typename NonlinearOptimizer<VectorType, MatrixType>::Constraints>
      inequality;
    if (inequality_matrix != nullptr)
      {
        inequality.emplace();
        inequality->values = *inequality_right_hand_side;
        Mat inequality_jacobian;
        internal::assert_nonlinear_optimizer_petsc(
          MatDuplicate(inequality_matrix->petsc_matrix(),
                       MAT_COPY_VALUES,
                       &inequality_jacobian));
        internal::reinit_linear_program_matrix(inequality->jacobian,
                                               inequality_jacobian,
                                               0);
        internal::assert_nonlinear_optimizer_petsc(
          MatDestroy(&inequality_jacobian));
        optimizer.inequality_constraints =
          [inequality_matrix, inequality_right_hand_side](const VectorType &x,
                                                          VectorType &values) {
            internal::assert_nonlinear_optimizer_petsc(
              MatMult(inequality_matrix->petsc_matrix(),
                      static_cast<const Vec &>(x),
                      values.petsc_vector()));
            internal::assert_nonlinear_optimizer_petsc(
              VecAXPY(values.petsc_vector(),
                      -1.,
                      static_cast<const Vec &>(*inequality_right_hand_side)));
          };
        optimizer.inequality_jacobian =
          [inequality_matrix](const VectorType &,
                              MatrixType &jacobian,
                              MatrixType &) {
            internal::assert_nonlinear_optimizer_petsc(
              MatCopy(inequality_matrix->petsc_matrix(),
                      jacobian.petsc_matrix(),
                      SAME_NONZERO_PATTERN));
          };
      }

    // Call the nonlinear optimizer's solve function with the constructed
    // callbacks and Hessian.
    return optimizer.solve(
      solution, hessian_wrapper, equality, inequality, bounds);
  }



  template <typename VectorType, typename MatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
  Tao LinearOptimizer<VectorType, MatrixType>::petsc_tao()
  {
    return optimizer.petsc_tao();
  }


} // namespace PETScWrappers
#  endif
#endif

DEAL_II_NAMESPACE_CLOSE

#endif
