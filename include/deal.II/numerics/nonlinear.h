// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_nonlinear_h
#define dealii_nonlinear_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi_stub.h>

#include <deal.II/lac/petsc_snes.h>
#include <deal.II/lac/vector.h>

#include <deal.II/sundials/kinsol.h>

#include <deal.II/trilinos/nox.h>


DEAL_II_NAMESPACE_OPEN

/**
 * This class offers a unified interface to several nonlinear solver
 * implementations like KINSOL, NOX, and SNES.
 *
 * KINSOL nonlinear solvers are part of the SUNDIALS package, while
 * NOX is part of Trilinos and SNES is part of PETSc, respectively.
 * If no solver is manually specified, this class will automaticlaly
 * choose one of the available solvers depending on the enabled
 * dependencies.
 *
 * By calling the @p solve function of this @p
 * NonlinearSolverSelector, it selects the @p solve function of that
 * @p Solver that was specified in the constructor of this class,
 * similar to the SolverSelector class.
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
 * nonlinear_solver.reinit_vector = [&](Vector<double> &x) {...};
 *
 * nonlinear_solver.residual =
 *           [&](const Vector<double> &evaluation_point,
 *               Vector<double> &      residual) {...};
 *
 * nonlinear_solver.setup_jacobian =
 *           [&](const Vector<double> &current_u) {...};
 *
 * nonlinear_solver.solve_with_jacobian =
 *           [&](const Vector<double> &rhs,
 *               Vector<double> &      dst,
 *               const double tolerance) {...};
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
       * Newton iteration with line search.
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
      /*
       * Use the PETSc SNES solver.
       */
      petsc_snes
    };

    /**
     * Initialization parameters for NonlinearSolverSelector.
     *
     * @param solver_type Nonlinear solver type.
     * @param strategy Method of solving the nonlinear problem.
     * @param maximum_non_linear_iterations Maximum number of nonlinear
     *   iterations.
     * @param function_tolerance Absolute stopping tolerance for the norm
     *   of the residual $F(u)$.
     * @param relative_tolerance Relative stopping tolerance.
     * @param step_tolerance Tolerance for minimum scaled step length
     * @param anderson_subspace_size Size of the Anderson acceleration
     *   subspace, use 0 to disable.
     */
    AdditionalData(const SolverType       &solver_type = automatic,
                   const SolutionStrategy &strategy    = linesearch,
                   const unsigned int      maximum_non_linear_iterations = 200,
                   const double            function_tolerance            = 1e-8,
                   const double            relative_tolerance            = 1e-5,
                   const double            step_tolerance                = 0.0,
                   const unsigned int      anderson_subspace_size        = 0);

    /**
     * The type of nonlinear solver to use. If the default 'automatic' is used,
     * it will choose the first available package in the order KINSOL, NOX, or
     * SNES.
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
     * If set to zero, default values will be used.
     */
    double function_tolerance;

    /**
     * Relative $l_2$ tolerance of the residual to be reached.
     *
     * @note Solver terminates successfully if either the function
     * tolerance or the relative tolerance has been reached.
     */
    const double relative_tolerance;

    /**
     * A scalar used as a stopping tolerance on the minimum
     * scaled step length.
     *
     * If set to zero, default values will be used.
     */
    double step_tolerance;

    /**
     * The size of the subspace used with Anderson acceleration
     * in conjunction with Picard iteration.
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
   * Constructor, selecting the solver and other parameters specified in
   * @p additional_data.
   *
   * @deprecated Use the other constructor with MPI_Comm instead.
   */
  DEAL_II_DEPRECATED
  NonlinearSolverSelector(const AdditionalData &additional_data);

  /**
   * Constructor.
   *
   * @param additional_data NonlinearSolverSelector configuration data
   * @param mpi_communicator MPI communicator used by the nonlinear solver.
   */
  NonlinearSolverSelector(const AdditionalData &additional_data,
                          const MPI_Comm       &mpi_communicator);

  /**
   * Select a new nonlinear solver. All solver names used in this class are
   * all lower case.
   */
  void
  select(const typename AdditionalData::SolverType &type);

  /**
   * Set the generic additional data for all nonlinear solvers.
   */
  void
  set_data(const AdditionalData &additional_data);

/**
 * Set the additional data for NOX. See TrilinosWrappers::NOXSolver
 * for more information.
 */
#ifdef DEAL_II_TRILINOS_WITH_NOX
  void
  set_data(
    const typename TrilinosWrappers::NOXSolver<VectorType>::AdditionalData
                                               &additional_data,
    const Teuchos::RCP<Teuchos::ParameterList> &parameters =
      Teuchos::rcp(new Teuchos::ParameterList));
#endif

/**
 * Set the additional data for KINSOL. See SUNDIALS::KINSOL for more
 * information.
 */
#ifdef DEAL_II_WITH_SUNDIALS
  void
  set_data(const typename SUNDIALS::KINSOL<VectorType>::AdditionalData
             &additional_data);
#endif

/**
 * Set the additional data for SNES. See PETScWrappers::NonlinearSolverData
 * for more information.
 */
#ifdef DEAL_II_WITH_PETSC
  void
  set_data(const typename PETScWrappers::NonlinearSolverData &additional_data);
#endif

  /**
   * Solve the nonlinear system.
   *
   * The content of `initial_guess_and_solution` is used as an initial guess
   * and the final solution is stored in the same vector.
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
   * @note This variable represents a
   * @ref GlossUserProvidedCallBack "user provided callback".
   * See there for a description of how to deal with errors and other
   * requirements and conventions. Some of the underlying packages
   * used by this class can deal with recoverable exceptions, whereas
   * others cannot. As a consequence, if a callback
   * throws an exception of type RecoverableUserCallbackError, then this
   * exception may or may not be treated like any other exception.
   */
  std::function<void(const VectorType &src, VectorType &dst)> residual;

  /**
   * A function object that users may supply and that is intended to
   * prepare the linear solver for subsequent calls to
   * solve_with_jacobian).
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
   * @param current_u Current value of $u$
   *
   * @note This variable represents a
   * @ref GlossUserProvidedCallBack "user provided callback".
   * See there for a description of how to deal with errors and other
   * requirements and conventions. Some of the underlying packages
   * used by this class can deal with recoverable exceptions, whereas
   * others cannot. As a consequence, if a callback
   * throws an exception of type RecoverableUserCallbackError, then this
   * exception may or may not be treated like any other exception.
   */
  std::function<void(const VectorType &current_u)> setup_jacobian;

  /**
   * A function object that users may supply and that is intended to solve
   * a linear system with the Jacobian matrix.
   *
   * @note PETSc SNES does not provide us with a linear tolerance to solve
   *       linear system with. We will provide a value of 1e-6 in that case.
   *
   * @param[in] rhs The system right hand side to solve for.
   * @param[out] dst The solution of $J^{-1} * \texttt{src}$.
   * @param[in] tolerance The tolerance with which to solve the linear system
   *   of equations.
   *
   * @note This variable represents a
   * @ref GlossUserProvidedCallBack "user provided callback".
   * See there for a description of how to deal with errors and other
   * requirements and conventions. Some of the underlying packages
   * used by this class can deal with recoverable exceptions, whereas
   * others cannot. As a consequence, if a callback
   * throws an exception of type RecoverableUserCallbackError, then this
   * exception may or may not be treated like any other exception.
   */
  std::function<
    void(const VectorType &rhs, VectorType &dst, const double tolerance)>
    solve_with_jacobian;

private:
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
#ifdef DEAL_II_TRILINOS_WITH_NOX
  typename TrilinosWrappers::NOXSolver<VectorType>::AdditionalData
                                       additional_data_nox;
  Teuchos::RCP<Teuchos::ParameterList> parameters_nox =
    Teuchos::rcp(new Teuchos::ParameterList);
#endif

/**
 * PETSc SNES configuration data
 */
#ifdef DEAL_II_WITH_PETSC
  typename PETScWrappers::NonlinearSolverData additional_data_petsc_snes;
#endif

  /**
   * Solve with PETSc SNES. Internal functions specialized for PETSc
   * Vectors. This is necessary to only instantiate the SNES solvers
   * with PETSc Vectors, as this is the only vector type currently
   * supported.
   */
  void
  solve_with_petsc(VectorType &initial_guess_and_solution);
};



template <typename VectorType>
void
NonlinearSolverSelector<VectorType>::set_data(
  const AdditionalData &additional_data)
{
  (void)additional_data;

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
#ifdef DEAL_II_TRILINOS_WITH_NOX
  parameters_nox->set("Nonlinear Solver", "Line Search Based");
  Teuchos::ParameterList &Line_Search = parameters_nox->sublist("Line Search");
  Line_Search.set("Method", "Full Step");

  additional_data_nox.max_iter = additional_data.maximum_non_linear_iterations;
  additional_data_nox.abs_tol  = additional_data.function_tolerance;
  additional_data_nox.rel_tol  = additional_data.relative_tolerance;
#endif

#ifdef DEAL_II_WITH_PETSC
  additional_data_petsc_snes.options_prefix = "";

  if (additional_data.anderson_subspace_size > 0 &&
      additional_data.strategy ==
        NonlinearSolverSelector<VectorType>::AdditionalData::picard)
    {
      additional_data_petsc_snes.snes_type = "anderson";
      // TODO additional_data.anderson_subspace_size;
    }
  else if (additional_data.strategy ==
           NonlinearSolverSelector<VectorType>::AdditionalData::linesearch)
    {
      additional_data_petsc_snes.snes_type            = "newtonls";
      additional_data_petsc_snes.snes_linesearch_type = "bt";
    }
  else if (additional_data.strategy ==
           NonlinearSolverSelector<VectorType>::AdditionalData::newton)
    {
      additional_data_petsc_snes.snes_linesearch_type = "newtonls";
      additional_data_petsc_snes.snes_linesearch_type = "basic";
    }
  else if (additional_data.strategy ==
           NonlinearSolverSelector<VectorType>::AdditionalData::picard)
    additional_data_petsc_snes.snes_type = "nrichardson";

  additional_data_petsc_snes.absolute_tolerance =
    additional_data.function_tolerance;
  additional_data_petsc_snes.relative_tolerance =
    additional_data.relative_tolerance;
  additional_data_petsc_snes.step_tolerance = additional_data.step_tolerance;
  additional_data_petsc_snes.maximum_non_linear_iterations =
    additional_data.maximum_non_linear_iterations;
  additional_data_petsc_snes.max_n_function_evaluations = -1;

#endif
}



template <typename VectorType>
NonlinearSolverSelector<VectorType>::NonlinearSolverSelector()
  : mpi_communicator(MPI_COMM_SELF)
{}



template <typename VectorType>
NonlinearSolverSelector<VectorType>::NonlinearSolverSelector(
  const AdditionalData &additional_data)
  : additional_data(additional_data)
  , mpi_communicator(MPI_COMM_SELF)
{
  set_data(additional_data);
}



template <typename VectorType>
NonlinearSolverSelector<VectorType>::NonlinearSolverSelector(
  const AdditionalData &additional_data,
  const MPI_Comm       &mpi_communicator)
  : additional_data(additional_data)
  , mpi_communicator(mpi_communicator)
{
  set_data(additional_data);
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
  const SolverType       &solver_type,
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



#ifdef DEAL_II_TRILINOS_WITH_NOX
template <typename VectorType>
void
NonlinearSolverSelector<VectorType>::set_data(
  const typename TrilinosWrappers::NOXSolver<VectorType>::AdditionalData
                                             &additional_data,
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
NonlinearSolverSelector<VectorType>::solve_with_petsc(
  VectorType & /*initial_guess_and_solution*/)
{
  AssertThrow(false,
              ExcMessage("PETSc SNES requires you to use PETSc vectors."));
}



#ifdef DEAL_II_WITH_PETSC
template <>
void
NonlinearSolverSelector<PETScWrappers::MPI::Vector>::solve_with_petsc(
  PETScWrappers::MPI::Vector &initial_guess_and_solution)
{
  PETScWrappers::NonlinearSolver<PETScWrappers::MPI::Vector> nonlinear_solver(
    additional_data_petsc_snes, mpi_communicator);

  nonlinear_solver.residual = residual;

  nonlinear_solver.setup_jacobian = setup_jacobian;

  nonlinear_solver.solve_with_jacobian =
    [&](const PETScWrappers::MPI::Vector &src,
        PETScWrappers::MPI::Vector       &dst) {
      // PETSc does not gives a tolerance, so we have to choose something
      // reasonable to provide to the user:
      const double tolerance = 1e-6;
      solve_with_jacobian(src, dst, tolerance);
    };

  nonlinear_solver.solve(initial_guess_and_solution);
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
#ifdef DEAL_II_WITH_PETSC
      additional_data.solver_type = AdditionalData::SolverType::petsc_snes;
#endif
#ifdef DEAL_II_TRILINOS_WITH_NOX
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

      nonlinear_solver.reinit_vector = reinit_vector;
      nonlinear_solver.residual      = residual;

      // We cannot simply set these two functions equal to each other
      // because they have a different number of inputs.
      nonlinear_solver.setup_jacobian = [&](const VectorType &current_u,
                                            const VectorType & /*current_f*/) {
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
#ifdef DEAL_II_TRILINOS_WITH_NOX

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
  else if (additional_data.solver_type ==
           AdditionalData::SolverType::petsc_snes)
    {
      // Forward to internal function specializations, which throws for
      // non-supported vector types:
      solve_with_petsc(initial_guess_and_solution);
    }
  else
    {
      const std::string solvers = ""
#ifdef DEAL_II_WITH_SUNDIALS
                                  "kinsol\n"
#endif
#ifdef DEAL_II_TRILINOS_WITH_NOX
                                  "NOX\n"
#endif
#ifdef DEAL_II_WITH_PETSC
                                  "SNES\n"
#endif
        ;

      AssertThrow(false,
                  ExcMessage(
                    "Invalid nonlinear solver specified. "
                    "The solvers available in your installation are:\n" +
                    solvers));
    }
}

DEAL_II_NAMESPACE_CLOSE

#endif
