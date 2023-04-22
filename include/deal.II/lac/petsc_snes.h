//-----------------------------------------------------------
//
//    Copyright (C) 2023 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//---------------------------------------------------------------

#ifndef dealii_petsc_snes_h
#define dealii_petsc_snes_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/parameter_handler.h>
#  include <deal.II/base/smartpointer.h>

#  include <deal.II/lac/petsc_matrix_base.h>
#  include <deal.II/lac/petsc_precondition.h>
#  include <deal.II/lac/petsc_vector_base.h>

#  include <petscsnes.h>

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
     * @note All parameters values specified here can be overriden by
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
     * PETSc solver type.
     */
    std::string snes_type;

    /**
     * Linesearch type.
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
   * class VectorType : public Subscriptor
   *    ...
   *    explicit VectorType(Vec);
   *    ...
   *    Vec & petsc_vector();
   *    ...
   * @endcode
   *
   * @code
   * class MatrixType : public Subscriptor
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
   * In alternative, users can also provide the implementation of the
   * Jacobian. This can be accomplished in two ways:
   *  - PETSc style using NonlinearSolver::jacobian
   *  - deal.II style using NonlinearSolver::setup_jacobian and
   * NonlinearSolver::solve_with_jacobian.
   *
   * In case both approaches are implemented, the deal.II style
   * will be used.
   *
   * The second approach is more in style with the deal.II philosophy
   * and it also allows command line customization, like for example,
   * @code
   * ./myApp -snes_type newtontr -ksp_type cg
   * @endcode
   * in case the user wants to change the default nonlinear solver to
   * a trust region solver and iterate on the matrix-free tangent system with
   * CG, still using NonlinearSolver::solve_with_jacobian as a preconditioner.
   *
   * The first approach has instead the advantage that only the matrix assembly
   * procedure has to be provided, thus allowing quicker implementations and
   * faster turnaround for experimenting with linear solver preconditioning
   * configurations via command line customizations, like for example,
   * @code
   * ./myApp -ksp_type cg -pc_type gamg
   * @endcode
   * See NonlinearSolver::set_matrix and NonlinearSolver::set_matrices for
   * additional details.
   *
   * In case the nonlinear equations are derived from energy minimization
   * arguments, it may be beneficial to perform linesearch or test trust-region
   * model reductions using the energy functional. In such cases, users can
   * set an optional NonlinearSolver::energy callback.
   *
   * @ingroup PETScWrappers
   */
  template <typename VectorType  = PETScWrappers::VectorBase,
            typename PMatrixType = PETScWrappers::MatrixBase,
            typename AMatrixType = PMatrixType>
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
                    const MPI_Comm &           mpi_comm = PETSC_COMM_WORLD);

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
     * Reset the solver. It does not change the customization.
     */
    void
    reinit();

    /**
     * Reset solver. Change customization according to @p data.
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
     * that PETSc will use.
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
     */
    std::function<int(const VectorType &x, VectorType &res)> residual;

    /**
     * Callback for the computation of the Jacobian
     * $\dfrac{\partial F}{\partial x}$.
     */
    std::function<int(const VectorType &x, AMatrixType &A, PMatrixType &P)>
      jacobian;

    /**
     * Callback for monitoring the solution process.
     *
     * This function is called by NonlinearSolver at the beginning
     * of each time step. Input arguments are the current step number
     * and the current value for ||F(x)||.
     */
    std::function<int(const VectorType & x,
                      const unsigned int step_number,
                      const real_type    f_norm)>
      monitor;

    /**
     * Callback to set up the Jacobian system.
     *
     * This callback gives full control to users to set up the tangent
     * operator $\dfrac{\partial F}{\partial x}$.
     *
     * Solvers must be provided via NonlinearSolver::solve_with_jacobian.
     */
    std::function<int(const VectorType &x)> setup_jacobian;

    /**
     * Callback for the solution of the tangent system set up with
     * NonlinearSolver::setup_jacobian.
     *
     * This is used as a preconditioner inside the Krylov process.
     */
    std::function<int(const VectorType &src, VectorType &dst)>
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
     */
    std::function<int(const VectorType &x, real_type &energy_value)> energy;

  private:
    /**
     * The PETSc SNES object.
     */
    SNES snes;

    SmartPointer<AMatrixType, NonlinearSolver> A;
    SmartPointer<PMatrixType, NonlinearSolver> P;

    /**
     * This flag is used to support versions of PETSc older than 3.13.
     */
    bool need_dummy_assemble;
  };

} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

#endif
