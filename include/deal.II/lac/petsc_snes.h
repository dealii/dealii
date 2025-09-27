// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_petsc_snes_h
#define dealii_petsc_snes_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/observer_pointer.h>
#  include <deal.II/base/parameter_handler.h>

#  include <deal.II/lac/petsc_matrix_base.h>
#  include <deal.II/lac/petsc_precondition.h>
#  include <deal.II/lac/petsc_vector_base.h>

#  include <petscsnes.h>

#  include <exception>

#  if defined(DEAL_II_HAVE_CXX20)
#    include <concepts>
#  endif


DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  /**
   * Additional parameters that can be passed to the NonlinearSolver class.
   */
  class NonlinearSolverData
  {
  public:
    /**
     * Type that holds real-valued numbers.
     *
     * Used to represent norms.
     */
    using real_type = PetscReal;

    /**
     * Initialization parameters for NonlinearSolverData.
     *
     * Running parameters:
     *
     * @param options_prefix The string indicating the options prefix for command line customization.
     * @param snes_type The string indicating the PETSc SNES solver type.
     * @param snes_linesearch_type The string indicating the PETSc linesearch type.
     * @param absolute_tolerance Absolute error tolerance.
     * @param relative_tolerance Relative error tolerance.
     * @param step_tolerance Step tolerance.
     * @param maximum_non_linear_iterations Maximum number of iterations allowed.
     * @param max_n_function_evaluations Maximum number of function evaluations allowed.
     *
     * @note All parameters values specified here can be overridden by
     * command line choices.
     *
     * @ingroup PETScWrappers
     */
    NonlinearSolverData(
      // Running parameters
      const std::string &options_prefix                = "",
      const std::string &snes_type                     = "",
      const std::string &snes_linesearch_type          = "",
      const real_type    absolute_tolerance            = 0,
      const real_type    relative_tolerance            = 0,
      const real_type    step_tolerance                = 0,
      const int          maximum_non_linear_iterations = -1,
      const int          max_n_function_evaluations    = -1)
      : options_prefix(options_prefix)
      , snes_type(snes_type)
      , snes_linesearch_type(snes_linesearch_type)
      , absolute_tolerance(absolute_tolerance)
      , relative_tolerance(relative_tolerance)
      , step_tolerance(step_tolerance)
      , maximum_non_linear_iterations(maximum_non_linear_iterations)
      , max_n_function_evaluations(max_n_function_evaluations)
    {}

    /**
     * Import parameter values.
     */
    void
    add_parameters(ParameterHandler &prm);

    /**
     * Options database prefix.
     */
    std::string options_prefix;

    /**
     * PETSc nonlinear solver type. Valid options include "newtonls"
     * (Newton with line search), "newtontr" (Newton with Trust
     * Region), "nrichardson" (Picard), and many more. See
     * https://petsc.org/release/manualpages/SNES/SNESType/ for more
     * information.
     */
    std::string snes_type;

    /**
     * Linesearch type. Valid options include "bt" (backtracking) and
     * "basic" (no line search). See
     * https://petsc.org/release/manualpages/SNES/SNESLineSearchType/
     * for more information.
     */
    std::string snes_linesearch_type;

    /**
     * Absolute error tolerance for function evaluation.
     *
     * @note Non-positive values indicate to use PETSc's default.
     */
    real_type absolute_tolerance;

    /**
     * Relative error tolerance for function evaluation.
     *
     * @note Non-positive values indicate to use PETSc's default.
     */
    real_type relative_tolerance;

    /**
     * Step tolerance for solution update.
     *
     * @note Non-positive values indicate to use PETSc's default.
     */
    real_type step_tolerance;

    /**
     * Maximum number of nonlinear iterations allowed.
     *
     * @note Negative values indicate to use PETSc's default.
     */
    int maximum_non_linear_iterations;

    /**
     * Maximum number of function evaluations allowed.
     *
     * @note Negative values indicate to use PETSc's default.
     */
    int max_n_function_evaluations;
  };

  /**
   * Interface to PETSc SNES solver for nonlinear equations.
   * The SNES solver is described in the
   * [PETSc manual](https://petsc.org/release/manual/snes/).
   *
   * This class solves the nonlinear system of algebraic
   * equations $F(x) = 0$.
   *
   * The interface to PETSc is realized by means of std::function callbacks
   * like in the TrilinosWrappers::NOXSolver and SUNDIALS::KINSOL classes.
   *
   * NonlinearSolver supports any vector and matrix type having constructors and
   * methods:
   *
   * @code
   * class VectorType : public EnableObserverPointer
   *    ...
   *    explicit VectorType(Vec);
   *    ...
   *    Vec & petsc_vector();
   *    ...
   * @endcode
   *
   * @code
   * class MatrixType : public EnableObserverPointer
   *    ...
   *    explicit MatrixType(Mat);
   *    ...
   *    Mat & petsc_matrix();
   *    ...
   * @endcode
   *
   * In particular, the supported types are the ones that can *wrap*
   * PETSc's native vector and matrix classes, that are able to modify
   * them in place, and that can return PETSc native types when requested.
   *
   * To use the solvers the user needs to provide the implementation of $F$ via
   * the NonlinearSolver::residual callback.
   *
   * The default linearization procedure of a solver instantiated with
   * this class consists in using Jacobian-Free-Newton-Krylov; the action of
   * tangent matrices inside a linear solver process are approximated via
   * matrix-free finite-differencing of the nonlinear residual equations.
   * For details, consult the [PETSc
   * manual](https://petsc.org/release/manual/snes/#sec-nlmatrixfree).
   *
   * Users can also provide the implementations of the Jacobian. This can be
   * accomplished in two ways:
   *  - PETSc style using NonlinearSolver::jacobian
   *  - deal.II style using NonlinearSolver::setup_jacobian and
   *    NonlinearSolver::solve_with_jacobian.
   * The preconditioning matrix can be specified using
   * NonlinearSolver::set_matrix(). In case both approaches are implemented, the
   * deal.II style will be used.
   *
   * NonlinearSolver::set_matrices() must be used in case the user wants to
   * provide the iteration matrix of the tangent system in the deal.II style
   * approach, thus replacing the matrix-free linearization.
   *
   * The correctness of the constructed Jacobians passed using
   * NonlinearSolver::set_matrix() can be checked using
   * @code
   * ./myApp -snes_test_jacobian
   * @endcode
   * See NonlinearSolver::set_matrix() and NonlinearSolver::set_matrices() for
   * additional details.
   *
   * The deal.II style approach still allows command line customization, like
   * for example,
   * @code
   * ./myApp -snes_type newtontr -ksp_type cg
   * @endcode
   * in case the user wants to change the default nonlinear solver to
   * a trust region solver and iterate on the matrix-free tangent system with
   * CG, still using NonlinearSolver::solve_with_jacobian as a preconditioner.
   *
   * The PETSc style approach has instead the advantage that only the matrix
   * assembly procedure has to be implemented, thus allowing quicker
   * implementations and faster turnaround for experimenting with linear solver
   * preconditioning configurations via command line customizations, like for
   * example,
   * @code
   * ./myApp -ksp_type cg -pc_type gamg
   * @endcode
   *
   * In case the nonlinear equations are derived from energy minimization
   * arguments, it may be beneficial to perform linesearch or test trust-region
   * model reductions using the energy functional. In such cases, users can
   * set an optional NonlinearSolver::energy callback.
   *
   * @dealiiConceptRequires{(concepts::is_dealii_petsc_vector_type<VectorType>
   * || std::constructible_from<VectorType,Vec>) &&
   *   (concepts::is_dealii_petsc_matrix_type<PMatrixType> ||
   *    std::constructible_from<PMatrixType,Mat>) &&
   *   (concepts::is_dealii_petsc_matrix_type<AMatrixType> ||
   *    std::constructible_from<AMatrixType, Mat>)}
   *
   * @ingroup PETScWrappers
   */
  template <typename VectorType  = PETScWrappers::VectorBase,
            typename PMatrixType = PETScWrappers::MatrixBase,
            typename AMatrixType = PMatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<PMatrixType> ||
               std::constructible_from<
                 PMatrixType,
                 Mat>)&&(concepts::is_dealii_petsc_matrix_type<AMatrixType> ||
                         std::constructible_from<AMatrixType, Mat>))
  class NonlinearSolver
  {
  public:
    /**
     * Type that holds real-valued numbers.
     *
     * Used to represent norms.
     */
    using real_type = PetscReal;

    /**
     * Constructor.
     */
    NonlinearSolver(const NonlinearSolverData &data     = NonlinearSolverData(),
                    const MPI_Comm             mpi_comm = PETSC_COMM_WORLD);

    /**
     * Destructor.
     */
    ~NonlinearSolver();

    /**
     * Conversion operator to gain access to the underlying PETSc type. If you
     * do this, you cut this class off some information it may need, so this
     * conversion operator should only be used if you know what you do.
     */
    operator SNES() const;

    /**
     * Return the PETSc SNES object.
     */
    SNES
    petsc_snes();

    /**
     * Return the underlying MPI communicator.
     */
    MPI_Comm
    get_mpi_communicator() const;

    /**
     * Reset the solver, it does not change the customization.
     */
    void
    reinit();

    /**
     * Reset solver.
     * Change customization according to @p data.
     */
    void
    reinit(const NonlinearSolverData &data);

    /**
     * Set the preconditioning matrix only.
     *
     * When used with NonlinearSolver::setup_jacobian and
     * NonlinearSolver::solve_with_jacobian, PETSc will approximate the
     * linear system matrix-vector product using an internal matrix-free
     * representation.
     *
     * When used with NonlinearSolver::jacobian PETSc will use the same matrix
     * for both preconditioning and matrix-vector products.
     */
    void
    set_matrix(PMatrixType &P);

    /**
     * Set both the linear system matrix and the preconditioning matrix
     * that PETSc will use (can be the same matrix). In this case, the
     * Jacobian-Free-Newton-Krylov approach will not be used.
     */
    void
    set_matrices(AMatrixType &A, PMatrixType &P);

    /**
     * Solve the nonlinear system of equations $F(x) = 0$.
     *
     * This function returns the number of iterations.
     * The vector @p x must contain the initial guess.
     * Upon returning, the @p x vector contains the solution.
     */
    unsigned int
    solve(VectorType &x);

    /**
     * Solve the nonlinear system of equations $F(x) = 0$.
     *
     * This function returns the number of iterations.
     * The vector @p x must contain the initial guess.
     * Upon returning, the @p x vector contains the solution.
     *
     * Here we also set the matrix to precondition the tangent system.
     */
    unsigned int
    solve(VectorType &x, PMatrixType &P);

    /**
     * Solve the nonlinear system of equations $F(x) = 0$.
     *
     * This function returns the number of iterations.
     * The vector @p x must contain the initial guess.
     * Upon returning, the @p x vector contains the solution.
     *
     * Here we also set the matrices to describe and precondition
     * the tangent system.
     */
    unsigned int
    solve(VectorType &x, AMatrixType &A, PMatrixType &P);

    /**
     * Callback for the computation of the nonlinear residual $F(x)$.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions.
     */
    std::function<void(const VectorType &x, VectorType &res)> residual;

    /**
     * Callback for the computation of the Jacobian
     * $\dfrac{\partial F}{\partial x}$.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions.
     */
    std::function<void(const VectorType &x, AMatrixType &A, PMatrixType &P)>
      jacobian;

    /**
     * Callback for monitoring the solution process.
     *
     * This function is called by NonlinearSolver at the beginning
     * of each step. Input arguments are the current step number
     * and the current value for $||F(x)||$.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions.
     */
    std::function<void(const VectorType  &x,
                       const unsigned int step_number,
                       const real_type    f_norm)>
      monitor;

    /**
     * Callback for the set up of the Jacobian system.
     *
     * This callback gives full control to users to set up the tangent
     * operator $\dfrac{\partial F}{\partial x}$.
     *
     * Solvers must be provided via NonlinearSolver::solve_with_jacobian.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions.
     */
    std::function<void(const VectorType &x)> setup_jacobian;

    /**
     * Callback for the solution of the tangent system set up with
     * NonlinearSolver::setup_jacobian.
     *
     * This is used as a preconditioner inside the Krylov process.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions.
     */
    std::function<void(const VectorType &src, VectorType &dst)>
      solve_with_jacobian;

    /**
     * Callback for the computation of the energy function.
     *
     * This is usually not needed, since by default SNES assumes that the
     * objective function to be minimized is $\frac{1}{2} || F(x) ||^2 $.
     *
     * However, if the nonlinear equations are derived from energy arguments, it
     * may be useful to use this callback to perform linesearch or to test for
     * the reduction in a trust region step.
     *
     * The value of the energy function must be returned in @p energy_value.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions.
     */
    std::function<void(const VectorType &x, real_type &energy_value)> energy;

  private:
    /**
     * The PETSc SNES object.
     */
    SNES snes;

    /**
     * Pointers to the internal PETSc matrix objects.
     */
    ObserverPointer<AMatrixType, NonlinearSolver> A;
    ObserverPointer<PMatrixType, NonlinearSolver> P;

    /**
     * This flag is used to support versions of PETSc older than 3.13.
     */
    bool need_dummy_assemble;

    /**
     * A pointer to any exception that may have been thrown in user-defined
     * call-backs and that we have to deal after the KINSOL function we call
     * has returned.
     */
    mutable std::exception_ptr pending_exception;
  };

} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

#endif
