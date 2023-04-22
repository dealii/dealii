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

#ifndef dealii_petsc_ts_h
#define dealii_petsc_ts_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/parameter_handler.h>
#  include <deal.II/base/smartpointer.h>

#  include <deal.II/lac/petsc_matrix_base.h>
#  include <deal.II/lac/petsc_precondition.h>
#  include <deal.II/lac/petsc_vector_base.h>

#  include <petscts.h>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  /**
   * Additional parameters that can be passed to the TimeStepper class.
   */
  class TimeStepperData
  {
  public:
    /**
     * Type that holds real-valued numbers.
     *
     * Used to represent time and norms tolerances.
     */
    using real_type = PetscReal;

    /**
     * Initialization parameters for TimeStepper.
     *
     * Running parameters:
     *
     * @param options_prefix The string indicating the options prefix for command line customization.
     * @param ts_type The string indicating the PETSc solver type.
     * @param initial_time Initial simulation time.
     * @param final_time Final simulation time.
     * @param initial_step_size Initial step size.
     * @param max_steps Maximum number of steps allowed.
     * @param match_step Whether or not to exactly stop at final time or step over it.
     *
     * Error parameters:
     *
     * @param ts_adapt_type The string indicating the PETSc time step adaptor type.
     * @param minimum_step_size Minimum step size allowed.
     * @param maximum_step_size Maximum step size allowed.
     * @param absolute_tolerance Absolute error tolerance.
     * @param relative_tolerance Relative error tolerance.
     * @param ignore_algebraic_lte Ignore algebraic terms for
     * error computations
     *
     * Note that one between `final_time` or `max_steps` must
     * be specified by the user, otherwise PETSc will complain.
     * Adaptive time stepping is disabled by default.
     * Negative values indicate using PETSc's default.
     *
     * @note All parameters values specified here can be overriden by
     * command line choices.
     *
     * @ingroup PETScWrappers
     */
    TimeStepperData(
      // Running parameters
      const std::string &options_prefix    = "",
      const std::string &ts_type           = "",
      const real_type    initial_time      = 0.0,
      const real_type    final_time        = 0.0,
      const real_type    initial_step_size = 0.0,
      const int          max_steps         = -1,
      const bool         match_step        = false,
      // Error parameters
      const std::string &ts_adapt_type        = "none",
      const real_type    minimum_step_size    = -1.0,
      const real_type    maximum_step_size    = -1.0,
      const real_type    absolute_tolerance   = -1.0,
      const real_type    relative_tolerance   = -1.0,
      const bool         ignore_algebraic_lte = true)
      : options_prefix(options_prefix)
      , ts_type(ts_type)
      , initial_time(initial_time)
      , final_time(final_time)
      , initial_step_size(initial_step_size)
      , max_steps(max_steps)
      , match_step(match_step)
      , ts_adapt_type(ts_adapt_type)
      , minimum_step_size(minimum_step_size)
      , maximum_step_size(maximum_step_size)
      , absolute_tolerance(absolute_tolerance)
      , relative_tolerance(relative_tolerance)
      , ignore_algebraic_lte(ignore_algebraic_lte)
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
    std::string ts_type;

    /**
     * Initial time for the DAE.
     */
    real_type initial_time;

    /**
     * Final time.
     */
    real_type final_time;

    /**
     * Initial step size.
     *
     * @note Non-positive values are ignored.
     */
    real_type initial_step_size;

    /**
     * Maximum number of steps to be taken.
     *
     * @note Negative values are ignored.
     */
    int max_steps;

    /**
     * Flag to indicate to stop exactly at the requested final time.
     */
    bool match_step;

    /**
     * PETSc time step adaptor type.
     */
    std::string ts_adapt_type;

    /**
     * Minimum allowed step size for adaptive time stepping.
     *
     * @note Non-positive values indicate to use PETSc's default.
     */
    real_type minimum_step_size;

    /**
     * Maximum allowed step size for adaptive time stepping.
     *
     * @note Non-positive values indicate to use PETSc's default.
     */
    real_type maximum_step_size;

    /**
     * Absolute error tolerance for adaptive time stepping.
     *
     * @note Negative values indicate to use PETSc's default.
     */
    real_type absolute_tolerance;

    /**
     * Relative error tolerance for adaptive time stepping.
     *
     * @note Negative values indicate to use PETSc's default.
     */
    real_type relative_tolerance;

    /**
     * Ignore algebraic terms for local truncation error.
     */
    bool ignore_algebraic_lte;
  };

  /**
   * Interface to the PETSc TS solver for Ordinary Differential Equations
   * and Differential-Algebraic Equations. The TS solver is described in the
   * [PETSc manual](https://petsc.org/release/manual/ts/).
   *
   * This class supports two kinds of formulations.
   * The explicit formulation:
   * \f[
   *   \begin{cases}
   *       \dot y = G(t,y)\, , \\
   *       y(t_0) = y_0\, , \\
   *   \end{cases}
   * \f]
   * and the implicit formulation:
   * \f[
   *   \begin{cases}
   *       F(t,y,\dot y) = 0\, , \\
   *       y(t_0) = y_0\, . \\
   *   \end{cases}
   * \f]
   *
   * The interface to PETSc is realized by means of std::function callbacks
   * like in the SUNDIALS::IDA and SUNDIALS::ARKode classes.
   *
   * TimeStepper supports any vector and matrix type having constructors and
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
   * To use explicit solvers (like for example explicit Runge-Kutta methods),
   * the user only needs to provide the implementation of $G$ via the
   * TimeStepper::explicit_function. For implicit solvers, users have also the
   * alternative of providing the $F$ function via
   * TimeStepper::implicit_function. IMEX methods are also supported by
   * providing both callbacks.
   *
   * The default linearization procedure of an implicit solver instantiated with
   * this class consists in using Jacobian-Free-Newton-Krylov; the action of
   * tangent matrices inside a linear solver process are approximated via
   * matrix-free finite-differencing of the nonlinear residual equations that
   * are ODE-solver specific. For details, consult the [PETSc
   * manual](https://petsc.org/release/manual/snes/#sec-nlmatrixfree).
   *
   * In alternative, users can also provide the implementations of the
   * *Jacobians*. This can be accomplished in two ways:
   *  - PETSc style using TimeStepper::implicit_jacobian
   *    and TimeStepper::explicit_jacobian.
   *  - deal.II style using TimeStepper::setup_jacobian and
   * TimeStepper::solve_with_jacobian
   *
   * In case both approaches are coded, the deal.II style
   * will be used.
   *
   * The second approach is more in style with the deal.II philosophy
   * and it also allows command line customization, like for example,
   * @code
   * ./myApp -snes_type newtontr -ksp_type cg
   * @endcode
   * in case the user wants to change the default nonlinear solver to
   * a trust region solver and iterate on the tangent system with CG,
   * still using TimeStepper::solve_with_jacobian as a preconditioner.
   *
   * The first approach has instead the advantage that only the matrix assembly
   * procedure has to be coded, thus allowing quicker implementations and
   * faster turnaround for experimenting with linear solver preconditioning
   * configurations via command line customizations, like for example,
   * @code
   * ./myApp -ksp_type cg -pc_type gamg
   * @endcode
   * See TimeStepper::set_matrix and TimeStepper::set_matrices for
   * additional details.
   *
   * @ingroup PETScWrappers
   */
  template <typename VectorType  = PETScWrappers::VectorBase,
            typename PMatrixType = PETScWrappers::MatrixBase,
            typename AMatrixType = PMatrixType>
  class TimeStepper
  {
  public:
    /**
     * Type that holds real-valued numbers.
     *
     * Used to represent time and norms tolerances.
     */
    using real_type = PetscReal;

    /**
     * Constructor.
     */
    TimeStepper(const TimeStepperData &data     = TimeStepperData(),
                const MPI_Comm &       mpi_comm = PETSC_COMM_WORLD);

    /**
     * Destructor.
     */
    ~TimeStepper();

    /**
     * Conversion operator to gain access to the underlying PETSc type. If you
     * do this, you cut this class off some information it may need, so this
     * conversion operator should only be used if you know what you do.
     */
    operator TS() const;

    /**
     * Return the PETSc TS object.
     */
    TS
    petsc_ts();

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
    reinit(const TimeStepperData &data);

    /**
     * Set the preconditioning matrix only.
     *
     * When used with TimeStepper::setup_jacobian and
     * TimeStepper::solve_with_jacobian, PETSc will approximate the linear
     * system matrix-vector product using an internal matrix-free
     * representation.
     *
     * When used with TimeStepper::implicit_jacobian or
     * TimeStepper::explicit_jacobian, PETSc will use the same matrix for both
     * preconditioning and matrix-vector products.
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
     * Return current time.
     */
    real_type
    get_time();

    /**
     * Return current time step.
     */
    real_type
    get_time_step();

    /**
     * Integrate the differential-algebraic equations starting from @p y.
     *
     * This function returns the final number of computed steps.
     * Upon returning, the @p y vector contains the solution of the DAE at
     * the end time.
     */
    unsigned int
    solve(VectorType &y);

    /**
     * Integrate the differential-algebraic equations starting from @p y.
     *
     * This function returns the final number of computed steps.
     * Upon returning, the @p y vector contains the solution of the DAE at
     * the end time.
     *
     * Here we also set the matrix to precondition the tangent system.
     */
    unsigned int
    solve(VectorType &y, PMatrixType &P);

    /**
     * Integrate the differential-algebraic equations starting from @p y.
     *
     * This function returns the final number of computed steps.
     * Upon returning, the @p y vector contains the solution of the DAE at
     * the end time.
     *
     * Here we also set the matrices to describe and precondition
     * the tangent system.
     */
    unsigned int
    solve(VectorType &y, AMatrixType &A, PMatrixType &P);

    /**
     * Callback for the computation of the implicit residual $F(t, y, \dot y)$.
     */
    std::function<int(const real_type   t,
                      const VectorType &y,
                      const VectorType &y_dot,
                      VectorType &      res)>
      implicit_function;

    /**
     * Callback for the computation of the implicit Jacobian
     * $\dfrac{\partial F}{\partial y} + \alpha \dfrac{\partial F}{\partial \dot
     * y}$.
     *
     * All implicit solvers implementations are recast to use the above
     * linearization. The $\alpha$ parameter is time-step and solver-type
     * specific.
     */
    std::function<int(const real_type   t,
                      const VectorType &y,
                      const VectorType &y_dot,
                      const real_type   alpha,
                      AMatrixType &     A,
                      PMatrixType &     P)>
      implicit_jacobian;

    /**
     * Callback for the computation of the explicit residual $G(t, y)$.
     */
    std::function<int(const real_type t, const VectorType &y, VectorType &res)>
      explicit_function;

    /**
     * Callback for the computation of the explicit Jacobian $\dfrac{\partial
     * G}{\partial y}$.
     */
    std::function<int(const real_type   t,
                      const VectorType &y,
                      AMatrixType &     A,
                      PMatrixType &     P)>
      explicit_jacobian;

    /**
     * Callback for monitoring the solution process.
     *
     * This function is called by TimeStepper at the beginning
     * of each time step.
     */
    std::function<int(const real_type    t,
                      const VectorType & y,
                      const unsigned int step_number)>
      monitor;

    /**
     * Callback for the set up of the Jacobian system.
     *
     * This callback gives full control to users to set up the linearized
     * equations
     * $\dfrac{\partial F}{\partial y} + \alpha \dfrac{\partial F}{\partial \dot
     * y}$.
     *
     * All implicit solvers implementations are recast to use the above
     * linearization. The $\alpha$ parameter is time-step and solver-type
     * specific.
     *
     * Solvers must be provided via TimeStepper::solve_with_jacobian.
     */
    std::function<int(const real_type   t,
                      const VectorType &y,
                      const VectorType &ydot,
                      const real_type   alpha)>
      setup_jacobian;

    /**
     * Callback for the solution of the tangent system set up with
     * TimeStepper::setup_jacobian.
     *
     * This is used as a preconditioner inside the Krylov process.
     */
    std::function<int(const VectorType &src, VectorType &dst)>
      solve_with_jacobian;

    /**
     * Callback to return an index set containing the algebraic components.
     *
     * Implementation of this function is optional. If your equation is also
     * algebraic (i.e., it contains algebraic constraints, or Lagrange
     * multipliers), you should implement this function in order to return only
     * these components of your system.
     */
    std::function<IndexSet()> algebraic_components;

  private:
    /**
     * The PETSc object.
     */
    TS ts;

    SmartPointer<AMatrixType, TimeStepper> A;
    SmartPointer<PMatrixType, TimeStepper> P;

    /**
     * This flag is set when changing the customization and used within solve.
     */
    bool need_dae_tolerances;

    /**
     * This flag is used to support versions of PETSc older than 3.13.
     */
    bool need_dummy_assemble;
  };

} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

#endif
