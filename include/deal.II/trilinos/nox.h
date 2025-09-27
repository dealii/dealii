// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_nox
#define dealii_trilinos_nox

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_NOX

#  include <deal.II/base/exceptions.h>

#  include <deal.II/lac/solver_control.h>

#  include <Teuchos_ParameterList.hpp>

#  include <exception>
#  include <functional>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  // Indicate that NOXSolver has not converged.
  DeclException0(ExcNOXNoConvergence);


  /**
   * Wrapper around the nonlinear solver from the NOX
   * package (https://docs.trilinos.org/dev/packages/nox/doc/html/index.html),
   * targeting deal.II data structures.
   *
   * The following code shows the steps how to use this class:
   * @code
   * // set configuration
   * TrilinosWrappers::NOXSolver<VectorType>::AdditionalData additional_data;
   *
   * // Define ParameterList object for more options
   * // These specifications are the default but we include them for
   * // clarification
   * const Teuchos::RCP<Teuchos::ParameterList> parameters =
   *   Teuchos::rcp(new Teuchos::ParameterList);
   *
   * // Specify nonlinear solver type
   * parameters->set("Nonlinear Solver","Line Search Based");
   *
   * // Specify method of line search
   * parameters->sublist("Line Search").set("Method","Full Step");
   *
   * // Specify direction
   * parameters->sublist("Direction").set("Method","Newton")
   *
   * // create nonlinear solver
   * TrilinosWrappers::NOXSolver<VectorType> solver(additional_data,parameters);
   *
   * // Set user functions to compute residual, to set up the Jacobian, and to
   * // apply the inverse of the Jacobian.
   * // Note that there are more functions that can be set.
   * solver.residual = [](const auto &u, auto &F) {...};
   * solver.setup_jacobian = [](const auto &u) {...};
   * solver.solve_with_jacobian =
   *   [](const auto &u, auto &F, const auto) {...};
   *
   * // solver nonlinear system with solution containing the initial guess and
   * // the final solution
   * solver.solve(solution);
   * @endcode
   *
   * The functions used in NOX are nearly identical to the functions in
   * SUNDIALS::KINSOL with a few exceptions (for example,
   * SUNDIALS::KINSOL requires a reinit() function where NOX does
   * not). So check the SUNDIALS::KINSOL documentation for more precise details
   * on how these functions are implemented.
   */
  template <typename VectorType>
  class NOXSolver
  {
  public:
    /**
     * Struct that helps to configure NOXSolver. More advanced
     * parameters are passed to the constructor NOXSolver
     * directly via a Teuchos::ParameterList.
     */
    struct AdditionalData
    {
    public:
      /**
       * Constructor.
       */
      AdditionalData(const unsigned int max_iter                       = 10,
                     const double       abs_tol                        = 1.e-20,
                     const double       rel_tol                        = 1.e-5,
                     const unsigned int threshold_nonlinear_iterations = 1,
                     const unsigned int threshold_n_linear_iterations  = 0,
                     const bool         reuse_solver                   = false);

      /**
       * Max number of nonlinear iterations.
       */
      unsigned int max_iter;

      /**
       * Absolute l2 tolerance of the residual to be reached.
       *
       * @note Solver terminates successfully if either the absolute or
       * the relative tolerance has been reached.
       */
      double abs_tol;

      /**
       * Relative l2 tolerance of the residual to be reached.
       *
       * @note Solver terminates successfully if either the absolute or
       * the relative tolerance has been reached.
       */
      double rel_tol;

      /**
       * Number of nonlinear iterations after which the preconditioner
       * should be updated.
       */
      unsigned int threshold_nonlinear_iterations;

      /**
       * A number that indicates how many iterations a linear solver
       * should at most perform before the preconditioner should
       * be updated. The use of this variable is predicated on the
       * idea that one can keep using a preconditioner built earlier
       * as long as it is a good preconditioner for the matrix currently
       * in use -- where "good" is defined as leading to a number of
       * iterations to solve linear systems less than the threshold
       * given by the current variable.
       *
       * This variable is only used if the
       * NOXSolver::solve_with_jacobian_and_track_n_linear_iterations
       * function object has been given a target (i.e., it is not empty).
       */
      unsigned int threshold_n_linear_iterations;

      /**
       * Reuse nonlinear solver the next time solve() is called. In particular,
       * this parameter allows to reuse the preconditioner from the last
       * solution step, enabling preconditioner lagging over multiple nonlinear
       * solutions.
       */
      bool reuse_solver;
    };

    /**
     * Constructor, with class parameters set by the AdditionalData object.
     *
     * @param additional_data NOX configuration data.
     * @param parameters More specific NOX solver configuration.
     *
     * If @p parameters is not filled, a Newton solver with a full step is used.
     * An overview of possible parameters is given at
     * https://docs.trilinos.org/dev/packages/nox/doc/html/parameters.html.
     */
    NOXSolver(AdditionalData                             &additional_data,
              const Teuchos::RCP<Teuchos::ParameterList> &parameters =
                Teuchos::rcp(new Teuchos::ParameterList));

    /**
     * Destructor.
     */
    ~NOXSolver();

    /**
     * Clear the internal state.
     */
    void
    clear();

    /**
     * Solve the nonlinear problem and return the number of nonlinear
     * iterations executed.
     */
    unsigned int
    solve(VectorType &solution);

    /**
     * A function object that users should supply and that is intended to
     * compute the residual $F(u)$.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. NOX can not deal
     * with "recoverable" errors for this callback, so if it
     * throws an exception of type RecoverableUserCallbackError, then this
     * exception is treated like any other exception.
     */
    std::function<void(const VectorType &u, VectorType &F)> residual;

    /**
     * A user function that sets up the Jacobian, based on the
     * current solution @p current_u.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. NOX can not deal
     * with "recoverable" errors for this callback, so if it
     * throws an exception of type RecoverableUserCallbackError, then this
     * exception is treated like any other exception.
     */
    std::function<void(const VectorType &current_u)> setup_jacobian;

    /**
     * A user function that sets up the preconditioner for inverting
     * the Jacobian, based on the current solution @p current_u.
     *
     * @note This function is optional and is used when setup_jacobian is
     * called and the preconditioner needs to be updated (see
     * update_preconditioner_predicate and
     * AdditionalData::threshold_nonlinear_iterations).
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. NOX can not deal
     * with "recoverable" errors for this callback, so if it
     * throws an exception of type RecoverableUserCallbackError, then this
     * exception is treated like any other exception.
     */
    std::function<void(const VectorType &current_u)> setup_preconditioner;

    /**
     * A user function that applies the Jacobian $\nabla_u F(u)$ to
     * @p x and writes the result in @p y. The Jacobian to be used
     * (i.e., more precisely: the linearization point $u$ above) is
     * the one computed when the `setup_jacobian` function was last called.
     *
     * @note This function is optional and is used in the case of certain
     * configurations. For instance, this function is required if the
     * polynomial line search (`NOX::LineSearch::Polynomial`) is
     * chosen, whereas for the full step case (`NOX::LineSearch::FullStep`)
     * it won't be called.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. NOX can not deal
     * with "recoverable" errors for this callback, so if it
     * throws an exception of type RecoverableUserCallbackError, then this
     * exception is treated like any other exception.
     */
    std::function<void(const VectorType &x, VectorType &y)> apply_jacobian;

    /**
     * A user function that applies the inverse of the Jacobian
     * $[\nabla_u F(u)]^{-1}$ to
     * @p y and writes the result in @p x. The parameter @p tolerance
     * specifies the error reduction if an iterative solver is used
     * in applying the inverse matrix. The Jacobian to be used
     * (i.e., more precisely: the linearization point $u$ above) is
     * the one computed when the `setup_jacobian` function was last called.
     *
     * @note This function is optional and is used in the case of certain
     * configurations.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. NOX can deal with "recoverable"
     * errors for this callback, if the NOX parameter
     * "Newton/Rescue Bad Newton Solve" is set to @p true (which is, in
     * fact, its default value). If this parameters is set to @p true,
     * then exceptions of type RecoverableUserCallbackError are eaten for
     * this callback and NOX can safely proceed with a recovery step.
     * Exceptions of other types are still treated as "irrecoverable".
     */
    std::function<
      void(const VectorType &y, VectorType &x, const double tolerance)>
      solve_with_jacobian;

    /**
     * A user function that applies the inverse of the Jacobian
     * $[\nabla_u F(u)]^{-1}$ to
     * @p y, writes the result in @p x and returns the number of
     * linear iterations the linear solver needed.
     * The parameter @p tolerance species the error reduction if an
     * iterative solver is used. The Jacobian to be used
     * (i.e., more precisely: the linearization point $u$ above) is
     * the one computed when the `setup_jacobian` function was last called.
     *
     * @note This function is used if `solve_with_jacobian` is not
     *   provided. Its return value is compared again
     *   AdditionalFlags::threshold_n_linear_iterations; if it is
     *   larger, the preconditioner will be built before the next
     *   linear system is solved. The use of this approach is predicated on the
     * idea that one can keep using a preconditioner built earlier
     * as long as it is a good preconditioner for the matrix currently
     * in use -- where "good" is defined as leading to a number of
     * iterations to solve linear systems less than the threshold
     * given by the current variable.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. NOX can deal with "recoverable"
     * errors for this callback, if the NOX parameter
     * "Newton/Rescue Bad Newton Solve" is set to @p true (which is, in
     * fact, its default value). If this parameters is set to @p true,
     * then exceptions of type RecoverableUserCallbackError are eaten for
     * this callback and NOX can safely proceed with a recovery step.
     * Exceptions of other types are still treated as "irrecoverable".
     */
    std::function<
      int(const VectorType &y, VectorType &x, const double tolerance)>
      solve_with_jacobian_and_track_n_linear_iterations;

    /**
     * A user function that allows to check convergence in addition to
     * ones checking the l2-norm and the number of iterations (see
     * AdditionalData). It is run after each nonlinear iteration.
     *
     * The input are the current iteration number @p i, the l2-norm
     * @p norm_f of the residual vector, the current solution @p current_u,
     * and the current residual vector @p f.
     *
     * @note This function is optional.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. NOX can not deal
     * with "recoverable" errors for this callback, so if it
     * throws an exception of type RecoverableUserCallbackError, then this
     * exception is treated like any other exception.
     */
    std::function<SolverControl::State(const unsigned int i,
                                       const double       norm_f,
                                       const VectorType  &current_u,
                                       const VectorType  &f)>
      check_iteration_status;

    /**
     * A user function that, in addition to
     * AdditionalData::threshold_nonlinear_iterations,
     * allows to force to update the preconditioner. A reason
     * for wanting to update the preconditioner is when the expected number
     * of linear iterations is exceeded.
     *
     * @note This function is optional. If no function is attached, this
     * means implicitly a return value of `false`.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. NOX can not deal
     * with "recoverable" errors for this callback, so if it
     * throws an exception of type RecoverableUserCallbackError, then this
     * exception is treated like any other exception.
     */
    std::function<bool()> update_preconditioner_predicate;

  private:
    /**
     * Additional data with basic settings.
     */
    AdditionalData additional_data;

    /**
     * Additional data with advanced settings. An overview of
     * possible parameters is given at
     * https://docs.trilinos.org/dev/packages/nox/doc/html/parameters.html.
     */
    const Teuchos::RCP<Teuchos::ParameterList> parameters;

    /**
     * A counter for the number of (accumulated) residual evaluations.
     */
    unsigned int n_residual_evaluations;

    /**
     * A counter for the number of (accumulated) Jacobi applications.
     */
    unsigned int n_jacobian_applications;

    /**
     * A counter for the number of (accumulated) nonlinear iterations.
     */
    unsigned int n_nonlinear_iterations;

    /**
     * The number of linear iterations of the last Jacobian solve.
     */
    unsigned int n_last_linear_iterations;

    /**
     * A pointer to any exception that may have been thrown in user-defined
     * call-backs and that we have to deal after the KINSOL function we call
     * has returned.
     */
    mutable std::exception_ptr pending_exception;
  };
} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif

#endif
