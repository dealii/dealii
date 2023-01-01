// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_trilinos_nox
#define dealii_trilinos_nox

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_NOX

#  include <deal.II/base/exceptions.h>

#  include <deal.II/lac/solver_control.h>

#  include <Teuchos_ParameterList.hpp>

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
   * // create nonlinear solver
   * TrilinosWrappers::NOXSolver<VectorType> solver(additional_data);
   *
   * // Set user functions to compute residual, to set up the Jacobian, and to
   * // apply the inverse of the Jacobian.
   * // Note that there are more functions that can be set.
   * solver.residual = [](const auto &src, auto &dst) {...};
   * solver.setup_jacobian = [](const auto &src) {...};
   * solver.solve_with_jacobian =
   *   [](const auto &src, auto &dst, const auto) {...};
   *
   * // solver nonlinear system with solution containing the initial guess and
   * // the final solution
   * solver.solve(solution);
   * @endcode
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
       * Max number of linear iterations after which the preconditioner
       * should be updated. This is only used if
       * solve_with_jacobian_and_track_n_linear_iterations has been given
       * a target (i.e., it is not empty).
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
     * Constructor.
     *
     * If @p parameters is not filled, a Newton solver with a full step is used.
     * An overview of possible parameters is given at
     * https://docs.trilinos.org/dev/packages/nox/doc/html/parameters.html.
     */
    NOXSolver(AdditionalData &                            additional_data,
              const Teuchos::RCP<Teuchos::ParameterList> &parameters =
                Teuchos::rcp(new Teuchos::ParameterList));

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
     * A user function that computes the residual @p f based on the
     * current solution @p x.
     *
     * @note This function should return 0 in the case of success.
     */
    std::function<int(const VectorType &x, VectorType &f)> residual;

    /**
     * A user function that sets up the Jacobian, based on the
     * current solution @p x.
     *
     * @note This function should return 0 in the case of success.
     */
    std::function<int(const VectorType &x)> setup_jacobian;

    /**
     * A user function that sets up the preconditioner for inverting
     * the Jacobian, based on the current solution @p x.
     *
     * @note This function is optional and is used when setup_jacobian is
     * called and the preconditioner needs to be updated (see
     * update_preconditioner_predicate and
     * AdditionalData::threshold_nonlinear_iterations).
     *
     * @note This function should return 0 in the case of success.
     */
    std::function<int(const VectorType &x)> setup_preconditioner;

    /**
     * A user function that applies the Jacobian to @p x and writes
     * the result in @p v.
     *
     * @note This function is optional and is used in the case of certain
     * configurations. For instance, this function is required if the
     * polynomial line search (`NOX::LineSearch::Polynomial`) is
     * chosen, whereas for the full step case (`NOX::LineSearch::FullStep`)
     * it won't be called.
     *
     * @note This function should return 0 in the case of success.
     */
    std::function<int(const VectorType &x, VectorType &v)> apply_jacobian;

    /**
     * A user function that applies the inverse of the Jacobian to
     * @p x and writes the result in @p x. The parameter @p tolerance
     * specifies the error reduction if an iterative solver is used.
     *
     * @note This function is optional and is used in the case of certain
     * configurations.
     *
     * @note This function should return 0 in the case of success.
     */
    std::function<
      int(const VectorType &f, VectorType &x, const double tolerance)>
      solve_with_jacobian;

    /**
     * A user function that applies the inverse of the Jacobian to
     * @p x, writes the result in @p x and returns the numer of
     * linear iterations the linear solver needed.
     * The parameter @p tolerance species the error reduction if a
     * interative solver is used.
     *
     * @note This function is optional and is used in the case of certain
     * configurations.
     */
    std::function<
      int(const VectorType &f, VectorType &x, const double tolerance)>
      solve_with_jacobian_and_track_n_linear_iterations;

    /**
     * A user function that allows to check convergence in addition to
     * ones checking the l2-norm and the number of iterations (see
     * AdditionalData). It is run after each nonlinear iteration.
     *
     * The input are the current iteration number @p i, the l2-norm
     * @p norm_f of the residual vector, the current solution @p x,
     * and the current residual vector @p f.
     *
     * @note This function is optional.
     */
    std::function<SolverControl::State(const unsigned int i,
                                       const double       norm_f,
                                       const VectorType & x,
                                       const VectorType & f)>
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
  };
} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
