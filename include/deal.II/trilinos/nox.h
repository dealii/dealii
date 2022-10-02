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

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/exceptions.h>

#  include <deal.II/lac/solver_control.h>

#  include <Teuchos_ParameterList.hpp>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  // Indicate that NOXSolver has not converged.
  DeclException0(ExcNOXNoConvergence);


  /**
   * Wrapper around the non-linear solver from the NOX
   * packge (https://docs.trilinos.org/dev/packages/nox/doc/html/index.html),
   * targeting deal.II data structures.
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
       * Max number of non-linear iterations.
       */
      unsigned int max_iter;

      /**
       * Absolute l2 tolerance to be reached.
       */
      double abs_tol;

      /**
       * Relative l2 tolerance to be reached.
       */
      double rel_tol;

      /**
       * Number of non-linear iterations after which the preconditioner
       * should be updated.
       */
      unsigned int threshold_nonlinear_iterations;

      /**
       * Max number of linear iterations after which the preconditioner
       * should be updated. This is only used if a lambda is attached to
       * solve_with_jacobian_and_track_n_linear_iterations.
       */
      unsigned int threshold_n_linear_iterations;

      /**
       * Reuse NOX solver instance in the next non-linear solution.
       */
      bool reuse_solver;
    };

    /**
     * Constructor.
     *
     * If @p parameters is not filled, a Newton solver with full step is used.
     * An overview of possible parameters is given at
     * https://docs.trilinos.org/dev/packages/nox/doc/html/parameters.html.
     */
    NOXSolver(AdditionalData &                            additional_data,
              const Teuchos::RCP<Teuchos::ParameterList> &parameters =
                Teuchos::rcp(new Teuchos::ParameterList));

    /**
     * Clear internal state.
     */
    void
    clear();

    /**
     * Solve non-linear problem and return number of iterations.
     */
    unsigned int
    solve(VectorType &solution);

    /**
     * User function that computes the residual.
     *
     * @note This function should return 0 in the case of success.
     */
    std::function<int(const VectorType &x, VectorType &f)> residual;

    /**
     * User function that sets up the Jacobian.
     *
     * @note This function should return 0 in the case of success.
     */
    std::function<int(const VectorType &x)> setup_jacobian;

    /**
     * User function that sets up the preconditioner for inverting
     * the Jacobian.
     *
     * @note The function is optional and is used when setup_jacobian is
     * called and the preconditioner needs to be updated (see
     * update_preconditioner_predicate and
     * AdditionalData::threshold_nonlinear_iterations).
     *
     * @note This function should return 0 in the case of success.
     */
    std::function<int(const VectorType &x)> setup_preconditioner;

    /**
     * User function that applies the Jacobian.
     *
     * @note The function is optional and is used in the case of certain
     * configurations. For instance, this function is required if the
     * polynomial line search (@p NOX::LineSearch::Polynomial) is
     * chosen, whereas for the full step case (@p NOX::LineSearch::FullStep)
     * it won't be called.
     *
     * @note This function should return 0 in the case of success.
     */
    std::function<int(const VectorType &x, VectorType &v)> apply_jacobian;

    /**
     * User function that applies the inverse of the Jacobian.
     *
     * @note The function is optional and is used in the case of certain
     * configurations.
     *
     * @note This function should return 0 in the case of success.
     */
    std::function<
      int(const VectorType &f, VectorType &x, const double tolerance)>
      solve_with_jacobian;

    /**
     * User function that applies the inverse of the Jacobian and
     * returns the numer of linear iterations the linear solver needed.
     *
     * @note This function should return -1 in the case of failure.
     */
    std::function<
      int(const VectorType &f, VectorType &x, const double tolerance)>
      solve_with_jacobian_and_track_n_linear_iterations;

    /**
     * User function that allows to check convergence in addition to
     * ones checking the l2-norm and the number of iterations (see
     * AdditionalData). It is run after each non-linear iteration.
     *
     * The input are the current iteration number @p i, the l2-norm
     * @p norm_f of the residual vector, the current solution @p x,
     * and the current residual vector @p f.
     *
     * @note The function is optional.
     */
    std::function<SolverControl::State(const unsigned int i,
                                       const double       norm_f,
                                       const VectorType & x,
                                       const VectorType & f)>
      check_iteration_status;

    /**
     * Function that allows to force to update the preconditioner in
     * addition to AdditionalData::threshold_nonlinear_iterations. A reason
     * for wanting to update the preconditioner is when the expected number
     * of linear iterations exceeds.
     *
     * @note The function is optional. If no function is attached, this
     * means implicitly a return value of false.
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
     * Counter for number of (accumulated) residual evaluations.
     */
    unsigned int n_residual_evaluations;

    /**
     * Counter for number of (accumulated) Jacobi applications.
     */
    unsigned int n_jacobian_applications;

    /**
     * Counter for number of (accumulated) non-linear iterations.
     */
    unsigned int n_nonlinear_iterations;

    /**
     * Number of linear iterations of the last Jacobian solve.
     */
    unsigned int n_last_linear_iterations;
  };
} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
