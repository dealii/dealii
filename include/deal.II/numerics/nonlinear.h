// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/sundials/kinsol.h>

#include <deal.II/trilinos/nox.h>


DEAL_II_NAMESPACE_OPEN

/**
 * Selects a nonlinear solver by choosing between KINSOL or NOX.
 * KINSOL and NOX are nonlinear solvers included in the SUNDIALS
 * package and the Trilinos package, respectively. If no solver type is
 * specified it will automaticlaly choose a solver based on what
 * deal.II was configured with, KINSOL having priority.
 *
 * By calling the @p solve function of this @p NonlinearSolverSelector, it selects the
 * @p solve function of that @p Solver that was specified in the constructor
 * of this class, similar to the SolverSelector class.
 *
 * <h3>Usage</h3>
 * An example of code one would run with this class is the
 * following:
 * @code
 * // Generate a @p NonlinearSolverSelector that uses @p KINSOL
 * NonlinearSolverSelector<Vector<double>>::AdditionalData additional_data;
 * additional_data.solver_type =
 *  NonlinearSolverSelector<Vector<double>>::AdditionalData::SolverType::kinsol;
 *
 * NonlinearSolverSelector<Vector<double>> nonlinear_solver(additional_data);
 *
 * // Functions that are required for solving a nonlinear problem,
 * // which are utilized in both @p KINSOL and @p NOX.
 * nonlinear_solver.reinit_vector = [&](Vector<double> &x) {...}
 *
 * nonlinear_solver.residual =
 *           [&](const Vector<double> &evaluation_point,
 *               Vector<double> &      residual) {...}
 *
 * nonlinear_solver.setup_jacobian =
 *           [&](const Vector<double> &current_u,
 *               const Vector<double> &current_f) {...}
 *
 * nonlinear_solver.solve_with_jacobian =
 *           [&](const Vector<double> &rhs,
 *               Vector<double> &      dst,
 *               const double tolerance) {...}
 *
 * // Calling the @p solve function with an initial guess.
 * nonlinear_solver.solve(current_solution);
 * @endcode
 */
template <typename VectorType = Vector<double>>
class NonlinearSolverSelector
{
public:
  /**
   * Additional data that will be sent through the NonlinearSolverSelector
   * class, then into the specific nonlinear solver classes themselves.
   * Many of these parameters can also be foud in the documentation for
   * @p KINSOL and @p NOX.
   */
  class AdditionalData
  {
  public:
    /**
     * NonlinearSolverSelector solution strategy. The solver types included in
     * this class both use a Newton-Krylov solver with a line search and even
     * Picard iterations (KINSOL).
     */
    enum SolutionStrategy
    {
      /**
       * Standard Newton iteration.
       */
      newton,
      /**
       * Newton iteration with linesearch.
       */
      linesearch,
      /**
       * Picard iteration.
       */
      picard,
    };

    enum SolverType
    {
      /**
       * Default parameter, will use whatever solver is available
       * with KINSOL as the priority.
       */
      automatic,
      /**
       * KINSOL nonlinear solver, part of the SUNDIALS package.
       */
      kinsol,
      /**
       * NOX nonlinear solver, part of the TRILINOS package.
       */
      nox,
    };

    /**
     * Initialization parameters for NonlinearSolverSelector.
     *
     * @param solver_type Nonlinear solver type, can be 'auto', 'kinsol', or 'nox'.
     * @param strategy Method of solving nonlinear problem, can be 'newton',
     * 'linesearch', or 'picard'.
     * @param maximum_non_linear_iterations Maximum number of nonlinear
     * iterations. This parameter is shared between KINSOL and NOX.
     * @param function_tolerance Function norm stopping tolerance.
     * @param relative_tolerance Relative function norm stopping tolerance.
     * @param step_tolerance Step tolerance for solution update.
     * @param anderson_subspace_size Anderson acceleration subspace size
     */
    AdditionalData(const SolverType &      solver_type = automatic,
                   const SolutionStrategy &strategy    = linesearch,
                   const unsigned int      maximum_non_linear_iterations = 200,
                   const double            function_tolerance            = 1e-8,
                   const double            relative_tolerance            = 1e-5,
                   const double            step_tolerance                = 0.0,
                   const unsigned int      anderson_subspace_size        = 0);

    /**
     * The type of nonlinear solver to use. The default value is set to 'auto',
     * which will use either KINSOL or NOX depending on which package is
     * installed, KINSOL having priority.
     */
    SolverType solver_type;

    /**
     * The solution strategy to use. For this class, you can choose from
     * SolutionStrategy::newton, SolutionStrategy::linesearch, or
     * SolutionStrategy::picard. More details on this can be found on the
     * @p KINSOL documentation.
     */
    SolutionStrategy strategy;

    /**
     * Maximum number of nonlinear iterations allowed.
     */
    unsigned int maximum_non_linear_iterations;

    /**
     * A scalar used as a stopping tolerance on the scaled
     * maximum norm of the system function $F(u)$ or $G(u)$.
     *
     * If set to zero, default values provided by KINSOL will be used.
     */
    double function_tolerance;

    /**
     * A scalar used as a stopping tolerance on the minimum
     * scaled step length.
     *
     * If set to zero, default values provided by KINSOL will be used.
     */
    double step_tolerance = 0.0;

    /**
     * Relative $l_2$ tolerance of the residual to be reached.
     *
     * @note Solver terminates successfully if either the absolute or
     * the relative tolerance has been reached.
     */
    const double relative_tolerance;

    /**
     * The size of the subspace used with Anderson acceleration
     * in conjunction with Picard or fixed-point iteration.
     *
     * If you set this to 0, no acceleration is used.
     */
    unsigned int anderson_subspace_size;
  };

  /**
   * Constructor, filling in default values
   */
  NonlinearSolverSelector();

  /**
   * Constructor, selecting the solver and other parametersspecified in
   * @p additional_data.
   */
  NonlinearSolverSelector(const AdditionalData &additional_data);

  /**
   * Constructor
   *
   * @param additional_data NonlinearSolverSelector configuration data
   * @param mpi_communicator MPI communicator over which logging operations are
   * computer.
   */
  NonlinearSolverSelector(const AdditionalData &additional_data,
                          const MPI_Comm &      mpi_communicator);

  /**
   * Select a new nonlinear solver. All solver names used in this class are
   * all lower case.
   */
  void
  select(const typename AdditionalData::SolverType &type);

/**
 * Set the additional data. For more information see the @p Solver class.
 */
#ifdef DEAL_II_WITH_TRILINOS
  void
  set_data(
    const typename TrilinosWrappers::NOXSolver<VectorType>::AdditionalData
      &                                         additional_data,
    const Teuchos::RCP<Teuchos::ParameterList> &parameters =
      Teuchos::rcp(new Teuchos::ParameterList));
#endif

/**
 * Set the additional data. For more information see the @p Solver class.
 */
#ifdef DEAL_II_WITH_SUNDIALS
  void
  set_data(const typename SUNDIALS::KINSOL<VectorType>::AdditionalData
             &additional_data);
#endif

  /**
   * Solve the nonlinear system. KINSOL uses the content of
   * `initial_guess_and_solution` as an initial guess, and
   * stores the final solution in the same vector.
   *
   * The functions herein are nearly identical in setup to what can be found
   * in the KINSOL and NOX documentation.
   */
  void
  solve(VectorType &initial_guess_and_solution);

  /**
   * A function object that users need to supply and that is intended to
   * reinitize the given vector to its correct size, and block structure (if
   * block vectors are used), along with any
   * other properties necessary.
   */
  std::function<void(VectorType &)> reinit_vector;

  /**
   * A function object that users should supply and that is intended to
   * compute the residual `dst = F(src)`.
   *
   * This function should return an int for either failure or success.
   */
  std::function<int(const VectorType &src, VectorType &dst)> residual;

  /**
   * A function object that users may supply and that is intended to
   * prepare the linear solver for subsequent calls to
   * solve_jacobian_system().
   *
   * The job of setup_jacobian() is to prepare the linear solver for
   * subsequent calls to solve_with_jacobian(), in the solution of linear
   * systems $Ax = b$. The exact nature of this system depends on the
   * SolutionStrategy that has been selected.
   *
   * In the cases strategy = SolutionStrategy::newton or
   * SolutionStrategy::linesearch, $A$ is the Jacobian $J = \partial
   * F/\partial u$. If strategy = SolutionStrategy::picard, $A$ is the
   * approximate Jacobian matrix $L$.
   *
   * The setup_jacobian() function may call a user-supplied function, or a
   * function within the linear solver module, to compute Jacobian-related
   * data that is required by the linear solver. It may also preprocess that
   * data as needed for solve_with_jacobian(), which may involve calling a
   * generic function (such as for LU factorization) or, more generally,
   * build preconditioners from the assembled Jacobian. In any case, the
   * data so generated may then be used whenever a linear system is solved.
   *
   * @param current_u Current value of $u$
   */
  std::function<int(const VectorType &current_u)> setup_jacobian;

  /**
   * A function object that users may supply and that is intended to solve
   * a linear system with the Jacobian matrix.
   *
   * Specific details on this function can be found in the KINSOL
   *
   * Arguments to the function are:
   *
   * @param[in] rhs The system right hand side to solve for.
   * @param[out] dst The solution of $J^{-1} * src$.
   * @param[in] tolerance The tolerance with which to solve the linear system
   *   of equations.
   */
  std::function<
    int(const VectorType &rhs, VectorType &dst, const double tolerance)>
    solve_with_jacobian;

protected:
  /**
   * NonlinearSolverSelector configuration data.
   */
  AdditionalData additional_data;

private:
  /**
   * The MPI communicator to be used by this solver, if any.
   */
  MPI_Comm mpi_communicator;

/**
 * KINSOL configuration data
 */
#ifdef DEAL_II_WITH_SUNDIALS
  typename SUNDIALS::KINSOL<VectorType>::AdditionalData additional_data_kinsol;
#endif

/**
 * NOX configuration data
 */
#ifdef DEAL_II_WITH_TRILINOS
  typename TrilinosWrappers::NOXSolver<VectorType>::AdditionalData
                                       additional_data_nox;
  Teuchos::RCP<Teuchos::ParameterList> parameters_nox =
    Teuchos::rcp(new Teuchos::ParameterList);
#endif

  /**
   * Data transfer function
   */
  void
  data_transfer(const AdditionalData &additional_data);
};

template <typename VectorType>
void
NonlinearSolverSelector<VectorType>::data_transfer(
  const AdditionalData &additional_data)
{
#ifdef DEAL_II_WITH_SUNDIALS
  // These if statements pass on the strategy to the other nonlinear solvers
  if (additional_data.strategy ==
      NonlinearSolverSelector<VectorType>::AdditionalData::linesearch)
    additional_data_kinsol.strategy =
      SUNDIALS::KINSOL<VectorType>::AdditionalData::linesearch;
  else if (additional_data.strategy ==
           NonlinearSolverSelector<VectorType>::AdditionalData::newton)
    additional_data_kinsol.strategy =
      SUNDIALS::KINSOL<VectorType>::AdditionalData::newton;
  else if (additional_data.strategy ==
           NonlinearSolverSelector<VectorType>::AdditionalData::picard)
    additional_data_kinsol.strategy =
      SUNDIALS::KINSOL<VectorType>::AdditionalData::picard;

  // Setting data points in the KINSOL class from the NonlinearSolverSelector
  // class
  additional_data_kinsol.maximum_non_linear_iterations =
    additional_data.maximum_non_linear_iterations;
  additional_data_kinsol.function_tolerance =
    additional_data.function_tolerance;
  additional_data_kinsol.step_tolerance = additional_data.step_tolerance;
  additional_data_kinsol.anderson_subspace_size =
    additional_data.anderson_subspace_size;
#endif

// Do the same thing we did above but with NOX
#ifdef DEAL_II_WITH_TRILINOS
  // Some default settings for parameters.
  parameters_nox->set("Nonlinear Solver", "Line Search Based");
  Teuchos::ParameterList &Line_Search = parameters_nox->sublist("Line Search");
  Line_Search.set("Method", "Full Step");

  additional_data_nox.max_iter = additional_data.maximum_non_linear_iterations;
  additional_data_nox.abs_tol  = additional_data.function_tolerance;
  additional_data_nox.rel_tol  = additional_data.relative_tolerance;
#endif
}

template <typename VectorType>
NonlinearSolverSelector<VectorType>::NonlinearSolverSelector() = default;

template <typename VectorType>
NonlinearSolverSelector<VectorType>::NonlinearSolverSelector(
  const AdditionalData &additional_data)
  : additional_data(additional_data)
{
  data_transfer(additional_data);
}

template <typename VectorType>
NonlinearSolverSelector<VectorType>::NonlinearSolverSelector(
  const AdditionalData &additional_data,
  const MPI_Comm &      mpi_communicator)
  : additional_data(additional_data)
  , mpi_communicator(mpi_communicator)
{
  data_transfer(additional_data);
}


template <typename VectorType>
void
NonlinearSolverSelector<VectorType>::select(
  const typename AdditionalData::SolverType &type)
{
  additional_data.solver_type = type;
}

template <typename VectorType>
NonlinearSolverSelector<VectorType>::AdditionalData::AdditionalData(
  const SolverType &      solver_type,
  const SolutionStrategy &strategy,
  const unsigned int      maximum_non_linear_iterations,
  const double            function_tolerance,
  const double            relative_tolerance,
  const double            step_tolerance,
  const unsigned int      anderson_subspace_size)
  : solver_type(solver_type)
  , strategy(strategy)
  , maximum_non_linear_iterations(maximum_non_linear_iterations)
  , function_tolerance(function_tolerance)
  , relative_tolerance(relative_tolerance)
  , step_tolerance(step_tolerance)
  , anderson_subspace_size(anderson_subspace_size)
{}

#ifdef DEAL_II_WITH_TRILINOS
template <typename VectorType>
void
NonlinearSolverSelector<VectorType>::set_data(
  const typename TrilinosWrappers::NOXSolver<VectorType>::AdditionalData
    &                                         additional_data,
  const Teuchos::RCP<Teuchos::ParameterList> &parameters)
{
  additional_data_nox = additional_data;
  parameters_nox      = parameters;
}
#endif

#ifdef DEAL_II_WITH_SUNDIALS
template <typename VectorType>
void
NonlinearSolverSelector<VectorType>::set_data(
  const typename SUNDIALS::KINSOL<VectorType>::AdditionalData &additional_data)
{
  additional_data_kinsol = additional_data;
}
#endif

template <typename VectorType>
void
NonlinearSolverSelector<VectorType>::solve(
  VectorType &initial_guess_and_solution)
{
  // The "auto" solver_type will default to kinsol, however if KINSOL is not
  // available then we will use NOX.
  if (additional_data.solver_type == AdditionalData::SolverType::automatic)
    {
#ifdef DEAL_II_WITH_TRILINOS
      additional_data.solver_type = AdditionalData::SolverType::nox;
#endif
#ifdef DEAL_II_WITH_SUNDIALS
      additional_data.solver_type = AdditionalData::SolverType::kinsol;
#endif

      // If "auto" is still the solver type we cannot solve the problem
      if (additional_data.solver_type == AdditionalData::SolverType::automatic)
        AssertThrow(false, ExcMessage("No valid solver type."));
    }

  if (additional_data.solver_type == AdditionalData::SolverType::kinsol)
    {
#ifdef DEAL_II_WITH_SUNDIALS
      SUNDIALS::KINSOL<VectorType> nonlinear_solver(additional_data_kinsol,
                                                    mpi_communicator);

      // We set the KINSOL reinit vector equal to the same function
      // defined for NonlinearSolverSelector.
      nonlinear_solver.reinit_vector = reinit_vector;

      nonlinear_solver.residual = residual;

      // We cannot simply set these two functions equal to each other
      // because they have a different number of inputs.
      nonlinear_solver.setup_jacobian = [&](const VectorType &current_u,
                                            const VectorType /*&current_f*/) {
        return NonlinearSolverSelector<VectorType>::setup_jacobian(current_u);
      };

      nonlinear_solver.solve_with_jacobian = solve_with_jacobian;

      nonlinear_solver.solve(initial_guess_and_solution);
#else
      AssertThrow(
        false, ExcMessage("You do not have SUNDIALS configured with deal.II!"));
#endif
    }
  else if (additional_data.solver_type == AdditionalData::SolverType::nox)
    {
#ifdef DEAL_II_WITH_TRILINOS
      TrilinosWrappers::NOXSolver<VectorType> nonlinear_solver(
        additional_data_nox, parameters_nox);

      // Do the same thing for NOX that we did with KINSOL.
      nonlinear_solver.residual = residual;

      // setup_jacobian for NOX has the same number of arguments for the same
      // function in NonlinearSolverSelector.
      nonlinear_solver.setup_jacobian = setup_jacobian;

      nonlinear_solver.solve_with_jacobian = solve_with_jacobian;

      nonlinear_solver.solve(initial_guess_and_solution);
#else
      AssertThrow(
        false, ExcMessage("You do not have Trilinos configured with deal.II"));
#endif
    }
  else
    {
      std::string solver1;
      std::string solver2;

#ifdef DEAL_II_WITH_SUNDIALS
      solver1 = "kinsol \n";
#endif
#ifdef DEAL_II_WITH_TRILINOS
      solver2 = "nox \n";
#endif

      DeclException2(InvalidNonlinearSolver,
                     std::string,
                     std::string,
                     "Invalid nonlinear solver specified, you may use:\n"
                       << arg1 << arg2);

      AssertThrow(false, InvalidNonlinearSolver(solver1, solver2));
    }
}


DEAL_II_NAMESPACE_CLOSE
