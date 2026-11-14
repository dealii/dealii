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

#ifndef dealii_petsc_linear_optimizer_h
#define dealii_petsc_linear_optimizer_h

#include <deal.II/base/config.h>

#include <deal.II/base/observer_pointer.h>

#if defined(DEAL_II_WITH_PETSC)
#  if DEAL_II_PETSC_VERSION_GTE(3, 13, 0)
#    include <deal.II/optimization/petsc_nonlinear_optimizer.h>
#  endif
#endif

DEAL_II_NAMESPACE_OPEN

#if defined(DEAL_II_WITH_PETSC)
#  if DEAL_II_PETSC_VERSION_GTE(3, 13, 0)


namespace PETScWrappers
{
  /**
   * Solver for linear programs using PETSc's TAO optimization framework. The
   * implementation of this class uses the NonlinearOptimizer class to solve
   * linear programs as a special case of nonlinear optimization problems. The
   * default solver, TAOALMM, uses an augmented Lagrangian multiplier method to
   * enforce the constraints, and TAOBNCG to solve its nested bound-constrained
   * subproblems.
   *
   * The NonlinearOptimizer class provides the necessary infrastructure to
   * handle the optimization process, including the management of constraints
   * and bounds, while the current class simplifies the interface for users
   * who are specifically interested in solving linear programs.
   *
   * This class solves linear programs of the form
   * \f[
   * \begin{aligned}
   *   \text{minimize} \quad & c^T x, \\
   *   \text{subject to} \quad & A_e x = b_e, \\
   *                            & A_i x - b_i \geq 0, \\
   *                            & \ell \leq x \leq u.
   * \end{aligned}
   * \f]
   * The equality and inequality constraints, as well as the variable bounds,
   * are optional.
   *
   * LinearOptimizer supports vector and matrix types that can wrap PETSc's
   * native `Vec` and `Mat` objects, modify them in place, and return the
   * corresponding PETSc object. This includes VectorBase and MatrixBase, as
   * well as PETScWrappers::MPI::Vector, PETScWrappers::MPI::SparseMatrix,
   * PETScWrappers::MPI::BlockVector, and
   * PETScWrappers::MPI::BlockSparseMatrix.
   *
   * VectorType and MatrixType must satisfy the same requirements as the
   * template arguments of NonlinearOptimizer.
   *
   * @dealiiConceptRequires{
   *   (concepts::is_dealii_petsc_vector_type<VectorType>
   *    || std::constructible_from<VectorType,Vec>) &&
   *   (concepts::is_dealii_petsc_matrix_type<MatrixType> ||
   *    std::constructible_from<MatrixType,Mat>)}
   *
   * <h3>Examples</h3>
   *
   * <h4>Unconstrained problems</h4>
   *
   * An unconstrained linear program has a finite minimum only when its
   * objective is zero. So this is perhaps not a particularly
   * relevant problem, but it illustrates the most basic use of
   * this class. So, to solve this kind of problem, use the following code:
   * @code
   * using VectorType = PETScWrappers::MPI::Vector;
   * using MatrixType = PETScWrappers::MPI::SparseMatrix;
   *
   * VectorType solution_vector; // Initialize to the initial guess.
   * VectorType objective;       // Initialize to zero for the current case.
   *
   * LinearOptimizer<VectorType, MatrixType> optimizer;
   * optimizer.solve(solution_vector, objective);
   * @endcode
   *
   *
   * <h4>Equality and inequality constraints</h4>
   *
   * To solve the linear program
   * \f[
   * \begin{aligned}
   *   \text{minimize} \quad & c^T x, \\
   *   \text{subject to} \quad & A_e x = b_e, \\
   *                            & A_i x \geq b_i,
   * \end{aligned}
   * \f]
   * where $A_e, A_i$ are matrices with as many columns as there are
   * variables, and $b_e, b_i$ are vectors with as many entries as the
   * corresponding matrices have rows, use:
   * @code
   * using VectorType = PETScWrappers::MPI::Vector;
   * using MatrixType = PETScWrappers::MPI::SparseMatrix;
   *
   * VectorType solution_vector; // Initialize to the initial guess.
   * VectorType objective;       // Initialize to c.
   * MatrixType equality_matrix; // Initialize to A_e.
   * VectorType equality_right_hand_side; // Initialize to b_e.
   * MatrixType inequality_matrix; // Initialize to A_i.
   * VectorType inequality_right_hand_side; // Initialize to b_i.
   *
   * LinearOptimizer<VectorType, MatrixType> optimizer;
   * optimizer.solve(solution_vector,
   *                 objective,
   *                 &equality_matrix,
   *                 &equality_right_hand_side,
   *                 &inequality_matrix,
   *                 &inequality_right_hand_side);
   * @endcode
   *
   *
   * <h4>Bound constraints</h4>
   *
   * To solve the bound-constraint linear program
   * \f[
   * \begin{aligned}
   *   \text{minimize} \quad & c^T x, \\
   *   \text{subject to} \quad & \ell \leq x \leq u,
   * \end{aligned}
   * \f]
   * where $\ell, u$ are vector describing element-wise bounds on the
   * solution vector $x$, use:
   * @code
   * using VectorType = PETScWrappers::MPI::Vector;
   * using MatrixType = PETScWrappers::MPI::SparseMatrix;
   * using Bounds =
   *   typename NonlinearOptimizer<VectorType, MatrixType>::Bounds;
   *
   * VectorType solution_vector; // Initialize to a feasible initial guess.
   * VectorType objective;       // Initialize to c.
   * Bounds bounds;              // Initialize bounds.lower and bounds.upper
   *                             // to the desired values. Set lower to
   *                             // PETSC_NINFINITY or upper to PETSC_INFINITY
   *                             // to disable a bound for an individual
   *                             // component.
   *
   * LinearOptimizer<VectorType, MatrixType> optimizer;
   * optimizer.solve(solution_vector,
   *                 objective,
   *                 nullptr,
   *                 nullptr,
   *                 nullptr,
   *                 nullptr,
   *                 bounds);
   * @endcode
   *
   * @ingroup PETScWrappers
   */
  template <typename VectorType = PETScWrappers::VectorBase,
            typename MatrixType = PETScWrappers::MatrixBase>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
  class LinearOptimizer
  {
  public:
    struct AdditionalData
      : public NonlinearOptimizer<VectorType, MatrixType>::AdditionalData
    {
      using Base =
        typename NonlinearOptimizer<VectorType, MatrixType>::AdditionalData;

      using ALMMSubsolver = typename Base::ALMMSubsolver;

      /**
       * Construct the solver parameters.
       *
       * @param options_prefix Prefix for PETSc options in the options database.
       * @param gradient_tolerance Absolute gradient tolerance. A non-positive
       * value selects PETSc's default.
       * @param constraint_tolerance Absolute constraint tolerance. A
       * non-positive value selects PETSc's default.
       * @param maximum_iterations Maximum number of optimization iterations. A
       * negative value selects PETSc's default.
       * @param almm_subsolver The solver for TAOALMM's nested
       * bound-constrained subproblem. Defaults to BNCG for linear programs.
       */
      AdditionalData(const std::string  &options_prefix       = "",
                     const PetscReal     gradient_tolerance   = 0.,
                     const PetscReal     constraint_tolerance = 0.,
                     const int           maximum_iterations   = -1,
                     const ALMMSubsolver almm_subsolver = ALMMSubsolver::bncg)
        : Base(options_prefix,
               gradient_tolerance,
               constraint_tolerance,
               maximum_iterations,
               almm_subsolver)
      {}
    };

    /**
     * Construct a linear-program solver.
     *
     * @param data Parameters controlling the optimization.
     * @param mpi_comm MPI communicator used by PETSc.
     */
    LinearOptimizer(const AdditionalData &data     = AdditionalData(),
                    const MPI_Comm        mpi_comm = PETSC_COMM_WORLD);

    /**
     * Copy constructor. Deleted to ensure that this object, which owns a PETSc
     * TAO object, cannot be copied.
     */
    LinearOptimizer(const LinearOptimizer &) = delete;

    /**
     * Copy assignment operator. Deleted to ensure that this object, which owns
     * a PETSc TAO object, cannot be copied.
     */
    LinearOptimizer &
    operator=(const LinearOptimizer &) = delete;

    /**
     * Solve a linear program. @p solution supplies the initial guess and is
     * overwritten with the solution. The constraints and bounds objects
     * are optional (see the class documentation for examples) but for
     * equality and inequality constraints, matrix and right-hand side must
     * either both be present or both be `nullptr`.
     *
     * @param solution Initial guess on entry and LP solution on return.
     * @param objective Coefficients of the linear objective.
     * @param equality_matrix Optional equality-constraint matrix.
     * @param equality_right_hand_side Optional equality-constraint right-hand
     * side.
     * @param inequality_matrix Optional inequality-constraint matrix.
     * @param inequality_right_hand_side Optional inequality-constraint
     * right-hand side.
     * @param bounds Optional lower and upper variable bounds.
     *
     * @return The number of optimization iterations performed.
     */
    unsigned int
    solve(VectorType                              &solution,
          const VectorType                        &objective,
          const ObserverPointer<MatrixType>       &equality_matrix = nullptr,
          const ObserverPointer<const VectorType> &equality_right_hand_side =
            nullptr,
          const ObserverPointer<MatrixType>       &inequality_matrix = nullptr,
          const ObserverPointer<const VectorType> &inequality_right_hand_side =
            nullptr,
          const std::optional<
            typename NonlinearOptimizer<VectorType, MatrixType>::Bounds>
            &bounds = std::nullopt);

    /**
     * Return the underlying PETSc TAO object.
     *
     * @return The PETSc TAO object.
     */
    Tao
    petsc_tao();

  private:
    /**
     * General nonlinear optimizer used to solve the linear program.
     */
    NonlinearOptimizer<VectorType, MatrixType> optimizer;
  };
} // namespace PETScWrappers
#  endif
#endif

DEAL_II_NAMESPACE_CLOSE

#endif
