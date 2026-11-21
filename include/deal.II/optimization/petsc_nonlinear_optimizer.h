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

#ifndef dealii_petsc_nonlinear_optimizer_h
#define dealii_petsc_nonlinear_optimizer_h

#include <deal.II/base/config.h>

#if defined(DEAL_II_WITH_PETSC)
#  if DEAL_II_PETSC_VERSION_GTE(3, 13, 0)
#    include <deal.II/base/mpi.h>
#    include <deal.II/base/observer_pointer.h>

#    include <deal.II/lac/petsc_matrix_base.h>
#    include <deal.II/lac/petsc_vector_base.h>
#    include <deal.II/lac/solver_control.h>

#    include <petsctao.h>

#    include <exception>
#    include <functional>
#    include <optional>
#    include <string>

#    if defined(DEAL_II_HAVE_CXX20)
#      include <concepts>
#    endif
#  endif
#endif

DEAL_II_NAMESPACE_OPEN

#if defined(DEAL_II_WITH_PETSC)
#  if DEAL_II_PETSC_VERSION_GTE(3, 13, 0)

namespace PETScWrappers
{
  /**
   * This class provides an interface to PETSc's nonlinear optimization
   * framework TAO for nonlinear objective functions, equality and inequality
   * constraints, and bound constraints. At present, the implementation uses
   * TAO's augmented Lagrangian multiplier method (TAOALMM), with TAOBQNLS as
   * its default nested subsolver.
   *
   * This class solves equality- and inequality-constrained optimization
   * problems of the form
   * \f[
   * \begin{aligned}
   *   \text{minimize} \quad & f(x), \\
   *   \text{subject to} \quad & g(x) = 0, \\
   *                            & h(x) \geq 0, \\
   *                            & \ell \leq x \leq u.
   * \end{aligned}
   * \f]
   * Here, $x\in\mathbb{R}^n$ is the optimization variable, $f$ is the
   * objective function, and $g$ and $h$ are vectors of equality and inequality
   * constraints, respectively. Equalities and inequalities are optional, as
   * are lower and upper bounds on the optimization; equals and greater-equals
   * signs above are understood to be componentwise.
   *
   * PETSc follows a caller-provided storage model: the user supplies vectors
   * and matrices that store the objective gradient and objective Hessian,
   * together with callbacks that update these objects at the current
   * optimization variable. Constraint residuals and Jacobians are owned by
   * the Constraints objects passed to solve(), and bounds are owned by the
   * Bounds object. PETSc invokes the callbacks as needed during the
   * optimization and reads their updated contents from the same supplied
   * storage.
   *
   * NonlinearOptimizer supports vector and matrix types that can wrap PETSc's
   * native `Vec` and `Mat` objects, modify them in place, and return the
   * corresponding PETSc object. This includes VectorBase and MatrixBase, as
   * well as PETScWrappers::MPI::Vector, PETScWrappers::MPI::SparseMatrix,
   * PETScWrappers::MPI::BlockVector, and
   * PETScWrappers::MPI::BlockSparseMatrix.
   *
   * @dealiiConceptRequires{
   *   (concepts::is_dealii_petsc_vector_type<VectorType>
   *    || std::constructible_from<VectorType,Vec>) &&
   *   (concepts::is_dealii_petsc_matrix_type<MatrixType> ||
   *    std::constructible_from<MatrixType,Mat>)}
   *
   *
   * <h3>Examples</h3>
   *
   * The following examples omit the bodies of the callbacks, which as
   * in other frameworks (for example the SUNDIALS::KINSOL or the
   * SUNDIALS::ARKode classes) are frequently implemented via lambda
   * functions. If the problem is unconstrained, callbacks for these
   * constraints do not have to be provided. However, the
   * objective, gradient, and Hessian callbacks must always be provided.
   *
   * <h4>Unconstrained problems</h4>
   *
   * To solve the unconstrained problem
   * \f[
   *   \text{minimize} \quad f(x),
   * \f]
   * use:
   * @code
   * using VectorType = PETScWrappers::MPI::Vector;
   * using MatrixType = PETScWrappers::MPI::SparseMatrix;
   *
   * // Set up the objects that will store the solution and Hessian.
   * VectorType solution_vector; // should be initialized to the initial guess
   * MatrixType hessian_matrix;
   *
   * // Now create the optimizer object and set the callbacks for the
   * // objective, gradient, and Hessian, and everything else as needed:
   * NonlinearOptimizer<VectorType, MatrixType> optimizer;
   * optimizer.objective = [&](const VectorType &x) -> PetscReal {
   *   // Compute the objective value at x.
   *   return value;
   * };
   * optimizer.gradient = [&](const VectorType &x,
   *                          VectorType       &gradient_vector) -> void {
   *   // Compute the objective gradient at x.
   * };
   * optimizer.hessian = [&](const VectorType &x,
   *                         MatrixType       &hessian_matrix,
   *                         MatrixType       &preconditioner) -> void {
   *   // Assemble the objective Hessian and preconditioner at x.
   * };
   * optimizer.solve(solution_vector, hessian_matrix);
   * @endcode
   *
   *
   * <h4>Equality and inequality constraints</h4>
   *
   * To solve the constrained problem
   * \f[
   * \begin{aligned}
   *   \text{minimize} \quad & f(x), \\
   *   \text{subject to} \quad & g(x) = 0, \\
   *                            & h(x) \geq 0,
   * \end{aligned}
   * \f]
   * represent each constraint set by a Constraints object, which owns its
   * residual vector and Jacobian:
   * @code
   * using VectorType = PETScWrappers::MPI::Vector;
   * using MatrixType = PETScWrappers::MPI::SparseMatrix;
   *
   * // Set up the objects that will store the solution and Hessian,
   * // along with vectors and Jacobians for the constraints:
   * VectorType solution_vector; // should be initialized to the initial guess
   * MatrixType hessian_matrix;
   *
   * // Now create the optimizer object and set the callbacks for the
   * // objective, gradient, and Hessian, and everything else as needed:
   * NonlinearOptimizer<VectorType, MatrixType>::Constraints equalities;
   * // Initialize equalities.values and equalities.jacobian.
   *
   * NonlinearOptimizer<VectorType, MatrixType>::Constraints inequalities;
   * // Initialize inequalities.values and inequalities.jacobian.
   *
   * NonlinearOptimizer<VectorType, MatrixType> optimizer;
   * optimizer.objective = [&](const VectorType &x) -> PetscReal {
   *   // Compute the objective value at x.
   *   return value;
   * };
   * optimizer.gradient = [&](const VectorType &x,
   *                          VectorType       &gradient_vector) -> void {
   *   // Compute the objective gradient at x.
   * };
   * optimizer.hessian = [&](const VectorType &x,
   *                         MatrixType       &hessian_matrix,
   *                         MatrixType       &preconditioner) -> void {
   *   // Assemble the objective Hessian and preconditioner at x.
   * };
   * optimizer.equality_constraints =
   *   [&](const VectorType &x, VectorType &values) -> void {
   *     // Compute g(x).
   *   };
   * optimizer.equality_jacobian =
   *   [&](const VectorType &x,
   *       MatrixType       &jacobian,
   *       MatrixType       &preconditioner) -> void {
   *     // Assemble the Jacobian of g at x.
   *   };
   * optimizer.inequality_constraints =
   *   [&](const VectorType &x, VectorType &values) -> void {
   *     // Compute h(x).
   *   };
   * optimizer.inequality_jacobian =
   *   [&](const VectorType &x,
   *       MatrixType       &jacobian,
   *       MatrixType       &preconditioner) -> void {
   *     // Assemble the Jacobian of h at x.
   *   };
   *
   * optimizer.solve(solution_vector,
   *                 hessian_matrix,
   *                 equalities,
   *                 inequalities);
   * @endcode
   *
   *
   * <h4>Bound constraints</h4>
   *
   * To solve the bound-constrained problem
   * \f[
   * \begin{aligned}
   *   \text{minimize} \quad & f(x), \\
   *   \text{subject to} \quad & \ell \leq x \leq u,
   * \end{aligned}
   * \f]
   * use:
   * @code
   * using VectorType = PETScWrappers::MPI::Vector;
   * using MatrixType = PETScWrappers::MPI::SparseMatrix;
   *
   * // Set up the objects that will store the solution and Hessian:
   * VectorType solution_vector; // should be initialized to a feasible initial
   * guess MatrixType hessian_matrix;
   *
   * NonlinearOptimizer<VectorType, MatrixType>::Bounds bounds;
   * // Initialize bounds.lower and bounds.upper.
   *
   * // Now create the optimizer object and set the callbacks for the
   * // objective, gradient, and Hessian, and everything else as needed:
   * NonlinearOptimizer<VectorType, MatrixType> optimizer;
   * optimizer.objective = [&](const VectorType &x) -> PetscReal {
   *   // Compute the objective value at x.
   *   return value;
   * };
   * optimizer.gradient = [&](const VectorType &x,
   *                          VectorType       &gradient_vector) -> void {
   *   // Compute the objective gradient at x.
   * };
   * optimizer.hessian = [&](const VectorType &x,
   *                         MatrixType       &hessian_matrix,
   *                         MatrixType       &preconditioner) -> void {
   *   // Assemble the objective Hessian and preconditioner at x.
   * };
   *
   * optimizer.solve(solution_vector,
   *                 hessian_matrix,
   *                 std::nullopt,
   *                 std::nullopt,
   *                 bounds);
   * @endcode
   *
   * @ingroup PETScWrappers
   */
  template <typename VectorType = PETScWrappers::VectorBase,
            typename MatrixType = PETScWrappers::MatrixBase>
#    if defined(DEAL_II_HAVE_CXX20) && \
      !defined(DEAL_II_DOXYGEN_DO_NOT_PARSE_REQUIRES_CLAUSES)
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<MatrixType> ||
               std::constructible_from<MatrixType, Mat>))
#    endif
  class NonlinearOptimizer
  {
  public:
    /**
     * Parameters controlling an optimization solve.
     */
    struct AdditionalData
    {
      /**
       * The method used for TAOALMM's nested bound-constrained subproblem.
       */
      enum class ALMMSubsolver
      {
        /**
         * PETSc's bounded quasi-Newton line-search method.
         */
        bqnls,

        /**
         * PETSc's bounded nonlinear conjugate-gradient method.
         */
        bncg
      };

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
       * bound-constrained subproblem.
       */
      AdditionalData(const std::string  &options_prefix       = "",
                     const PetscReal     gradient_tolerance   = 0.,
                     const PetscReal     constraint_tolerance = 0.,
                     const int           maximum_iterations   = -1,
                     const ALMMSubsolver almm_subsolver = ALMMSubsolver::bqnls)
        : options_prefix(options_prefix)
        , gradient_tolerance(gradient_tolerance)
        , constraint_tolerance(constraint_tolerance)
        , maximum_iterations(maximum_iterations)
        , almm_subsolver(almm_subsolver)
      {}

      /**
       * Prefix for PETSc options in the options database.
       */
      std::string options_prefix;

      /**
       * Absolute gradient tolerance.
       */
      PetscReal gradient_tolerance;

      /**
       * Absolute constraint tolerance.
       */
      PetscReal constraint_tolerance;

      /**
       * Maximum number of optimization iterations.
       */
      int maximum_iterations;

      /**
       * The solver for TAOALMM's nested bound-constrained subproblem.
       */
      ALMMSubsolver almm_subsolver;
    };

    /**
     * A constraint residual vector and its Jacobian matrix. This object owns
     * the storage passed to PETSc and must remain alive for the duration of
     * solve().
     */
    struct Constraints
    {
      /**
       * Storage for the constraint residual.
       */
      VectorType values;

      /**
       * Storage for the constraint Jacobian.
       */
      MatrixType jacobian;
    };

    /**
     * Lower and upper bounds for the optimization variables. Set an entry of
     * lower to `PETSC_NINFINITY` or an entry of upper to `PETSC_INFINITY` to
     * disable the corresponding bound.
     */
    struct Bounds
    {
      /**
       * Lower bounds for the optimization variables.
       */
      VectorType lower;

      /**
       * Upper bounds for the optimization variables.
       */
      VectorType upper;
    };

    /**
     * Construct a nonlinear optimizer.
     *
     * @param data Parameters controlling the optimization.
     * @param mpi_comm MPI communicator used by PETSc.
     */
    NonlinearOptimizer(const AdditionalData &data     = AdditionalData(),
                       const MPI_Comm        mpi_comm = PETSC_COMM_WORLD);

    /**
     * Copy constructor. Deleted to ensure that this object, which owns a PETSc
     * TAO object, cannot be copied.
     */
    NonlinearOptimizer(const NonlinearOptimizer &) = delete;

    /**
     * Copy assignment operator. Deleted to ensure that this object, which owns
     * a PETSc TAO object, cannot be copied.
     */
    NonlinearOptimizer &
    operator=(const NonlinearOptimizer &) = delete;

    /**
     * Destructor.
     */
    ~NonlinearOptimizer();

    /**
     * Solve the optimization problem. @p solution supplies the initial guess
     * and is replaced on exit from this function by the computed solution.
     * When @p bounds is supplied, each entry of the initial guess must satisfy
     * its lower and upper bound.
     *
     * @param solution Initial guess on entry and optimization solution on
     * return. When @p bounds is supplied, the initial guess must satisfy
     * those bounds.
     * @param hessian Storage for the Hessian of the objective.
     * @param equality_constraints Optional equality constraint residual and
     * Jacobian.
     * @param inequality_constraints Optional inequality constraint residual and
     * Jacobian.
     * @param bounds Optional lower and upper variable bounds.
     *
     * @return The number of optimization iterations performed.
     *
     * @throws SolverControl::NoConvergence if TAO terminates without
     * convergence within the allowed number of iterations. The exception
     * object contains the number of iterations performed and the final
     * residual norm (i.e., the norm of the gradient of the objective
     * function).
     */
    unsigned int
    solve(
      VectorType                       &solution,
      MatrixType                       &hessian,
      const std::optional<Constraints> &equality_constraints   = std::nullopt,
      const std::optional<Constraints> &inequality_constraints = std::nullopt,
      const std::optional<Bounds>      &bounds                 = std::nullopt);

    /**
     * Return the underlying PETSc TAO object. This can be used to set PETSc
     * options from the options database, for example, or to query the
     * optimization status after solve() returns.
     *
     * @return The PETSc TAO object.
     */
    Tao
    petsc_tao();

    /**
     * Callback that evaluates the objective function.
     *
     * @param x Current optimization variable.
     *
     * @return The objective value.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions.
     */
    std::function<PetscReal(const VectorType &)> objective;

    /**
     * Callback that evaluates the objective gradient.
     *
     * @param x Current optimization variable.
     * @param gradient Computed objective gradient.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions.
     */
    std::function<void(const VectorType &, VectorType &)> gradient;

    /**
     * Callback that assembles the Hessian of the objective function.
     *
     * @param x Current optimization variable.
     * @param hessian Hessian matrix to assemble.
     * @param preconditioner Preconditioning matrix to assemble.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions.
     */
    std::function<void(const VectorType &, MatrixType &, MatrixType &)> hessian;

    /**
     * Callback that evaluates the equality constraints.
     *
     * @param x Current optimization variable.
     * @param values Computed equality-constraint residual.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions.
     */
    std::function<void(const VectorType &, VectorType &)> equality_constraints;

    /**
     * Callback that assembles the equality-constraint Jacobian.
     *
     * @param x Current optimization variable.
     * @param jacobian Equality-constraint Jacobian to assemble.
     * @param preconditioner Preconditioning matrix to assemble.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions.
     */
    std::function<void(const VectorType &, MatrixType &, MatrixType &)>
      equality_jacobian;

    /**
     * Callback that evaluates the inequality constraints. PETSc uses the
     * convention $h(x) \geq 0$.
     *
     * @param x Current optimization variable.
     * @param values Computed inequality-constraint residual.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions.
     */
    std::function<void(const VectorType &, VectorType &)>
      inequality_constraints;

    /**
     * Callback that assembles the inequality-constraint Jacobian.
     *
     * @param x Current optimization variable.
     * @param jacobian Inequality-constraint Jacobian to assemble.
     * @param preconditioner Preconditioning matrix to assemble.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions.
     */
    std::function<void(const VectorType &, MatrixType &, MatrixType &)>
      inequality_jacobian;

  private:
    /**
     * Create and configure the PETSc TAO object.
     */
    void
    initialize_tao();

    /**
     * Forward PETSc's objective request to the user callback.
     *
     * @param tao PETSc TAO object.
     * @param x Current optimization variable.
     * @param value Storage for the objective value.
     * @param context Pointer to this optimizer.
     */
    static PetscErrorCode
    objective_callback(Tao tao, Vec x, PetscReal *value, void *context);

    /**
     * Forward PETSc's gradient request to the user callback.
     *
     * @param tao PETSc TAO object.
     * @param x Current optimization variable.
     * @param gradient Storage for the objective gradient.
     * @param context Pointer to this optimizer.
     */
    static PetscErrorCode
    gradient_callback(Tao tao, Vec x, Vec gradient, void *context);

    /**
     * Forward PETSc's Hessian request to the user callback.
     *
     * @param tao PETSc TAO object.
     * @param x Current optimization variable.
     * @param hessian Storage for the Hessian.
     * @param preconditioner Storage for the preconditioning matrix.
     * @param context Pointer to this optimizer.
     */
    static PetscErrorCode
    hessian_callback(Tao   tao,
                     Vec   x,
                     Mat   hessian,
                     Mat   preconditioner,
                     void *context);

    /**
     * Forward PETSc's equality-constraint request to the user callback.
     *
     * @param tao PETSc TAO object.
     * @param x Current optimization variable.
     * @param values Storage for the equality-constraint residual.
     * @param context Pointer to this optimizer.
     */
    static PetscErrorCode
    equality_constraints_callback(Tao tao, Vec x, Vec values, void *context);

    /**
     * Forward PETSc's equality-Jacobian request to the user callback.
     *
     * @param tao PETSc TAO object.
     * @param x Current optimization variable.
     * @param jacobian Storage for the equality-constraint Jacobian.
     * @param preconditioner Storage for the preconditioning matrix.
     * @param context Pointer to this optimizer.
     */
    static PetscErrorCode
    equality_jacobian_callback(Tao   tao,
                               Vec   x,
                               Mat   jacobian,
                               Mat   preconditioner,
                               void *context);

    /**
     * Forward PETSc's inequality-constraint request to the user callback.
     *
     * @param tao PETSc TAO object.
     * @param x Current optimization variable.
     * @param values Storage for the inequality-constraint residual.
     * @param context Pointer to this optimizer.
     */
    static PetscErrorCode
    inequality_constraints_callback(Tao tao, Vec x, Vec values, void *context);

    /**
     * Forward PETSc's inequality-Jacobian request to the user callback.
     *
     * @param tao PETSc TAO object.
     * @param x Current optimization variable.
     * @param jacobian Storage for the inequality-constraint Jacobian.
     * @param preconditioner Storage for the preconditioning matrix.
     * @param context Pointer to this optimizer.
     */
    static PetscErrorCode
    inequality_jacobian_callback(Tao   tao,
                                 Vec   x,
                                 Mat   jacobian,
                                 Mat   preconditioner,
                                 void *context);

    /**
     * PETSc TAO object owned by this optimizer.
     */
    Tao tao;

    /**
     * Parameters used to configure TAO.
     */
    const AdditionalData additional_data;

    /**
     * Communicator used to create TAO.
     */
    const MPI_Comm mpi_comm;

    /**
     * Whether a previous solve has configured TAO's constraints.
     */
    bool constraints_are_initialized;

    /**
     * Constraint sets supplied to the previous solve.
     */
    bool previous_has_equality_constraints;
    bool previous_has_inequality_constraints;

    /**
     * Exception captured in a user callback and rethrown after TAO returns.
     */
    std::exception_ptr pending_exception;
  };
} // namespace PETScWrappers
#  endif
#endif

DEAL_II_NAMESPACE_CLOSE

#endif
