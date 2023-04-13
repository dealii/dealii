//-----------------------------------------------------------
//
//    Copyright (C) 2022 by the deal.II authors
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
//
// Author: Stefano Zampini, King Abdullah University of Science and Technology.

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
     * @param opts_prefix The string indicating the options prefix for command line customization.
     * @param tstype The string indicating the PETSc solver type.
     * @param initial_time Initial simulation time.
     * @param final_time Final simulation time.
     * @param initial_step_size Initial step size.
     * @param max_steps Maximum number of steps allowed.
     * @param match_step Whether or not to exactly stop at final time or step over it.
     *
     * Error parameters:
     *
     * @param tsadapttype The string indicating the PETSc time step adaptor type.
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
     * Note that all parameters values specified here can be overriden by
     * command line choices.
     */
    TimeStepperData(
      // Running parameters
      const std::string &opts_prefix       = "",
      const std::string &tstype            = "",
      const real_type    initial_time      = 0.0,
      const real_type    final_time        = 0.0,
      const real_type    initial_step_size = 0.0,
      const int          max_steps         = -1,
      const bool         match_step        = false,
      // Error parameters
      const std::string &tsadapttype          = "none",
      const real_type    minimum_step_size    = -1.0,
      const real_type    maximum_step_size    = -1.0,
      const real_type    absolute_tolerance   = -1.0,
      const real_type    relative_tolerance   = -1.0,
      const bool         ignore_algebraic_lte = true)
      : opts_prefix(opts_prefix)
      , tstype(tstype)
      , initial_time(initial_time)
      , final_time(final_time)
      , initial_step_size(initial_step_size)
      , max_steps(max_steps)
      , match_step(match_step)
      , tsadapttype(tsadapttype)
      , minimum_step_size(minimum_step_size)
      , maximum_step_size(maximum_step_size)
      , absolute_tolerance(absolute_tolerance)
      , relative_tolerance(relative_tolerance)
      , ignore_algebraic_lte(ignore_algebraic_lte)
    {}

    void
    add_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Running parameters");
      prm.add_parameter(
        "options prefix",
        opts_prefix,
        "The string indicating the options prefix for command line customization.");
      prm.add_parameter("solver type",
                        tstype,
                        "The string indicating the PETSc TS type.");
      prm.add_parameter("initial time",
                        initial_time,
                        "The value for the initial time.");
      prm.add_parameter("final time",
                        final_time,
                        "The value for the final time.");
      prm.add_parameter("initial step size",
                        initial_step_size,
                        "The value for the initial time step.");
      prm.add_parameter("maximum number of steps",
                        max_steps,
                        "Maximum number of time steps allowed.");
      prm.add_parameter(
        "match final time",
        match_step,
        "Whether or not to exactly stop at final time or step over it.");
      prm.leave_subsection();

      prm.enter_subsection("Error control");
      prm.add_parameter("adaptor type",
                        tsadapttype,
                        "The string for the TSAdapt type.");
      prm.add_parameter("minimum step size",
                        minimum_step_size,
                        "Minimum time step size allowed.");
      prm.add_parameter("maximum step size",
                        maximum_step_size,
                        "Maximum time step size allowed.");
      prm.add_parameter("absolute error tolerance",
                        absolute_tolerance,
                        "Absolute error tolerance.");
      prm.add_parameter("relative error tolerance",
                        relative_tolerance,
                        "Absolute error tolerance.");
      prm.add_parameter(
        "ignore algebraic lte",
        ignore_algebraic_lte,
        "Indicate whether or not to suppress algebraic variables "
        "in the local truncation error test.");
      prm.leave_subsection();
    }

    /**
     * Options database prefix.
     */
    std::string opts_prefix;

    /**
     * PETSc solver type.
     */
    std::string tstype;

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
     * Non-positive values are ignored.
     */
    real_type initial_step_size;

    /**
     * Maximum number of steps to be taken.
     * Negative values are ignored.
     */
    int max_steps;

    /**
     * Match final time requested?
     */
    bool match_step;

    /**
     * PETSc time step adaptor type.
     */
    std::string tsadapttype;

    /**
     * Minimum allowed step size for adaptive time stepping.
     * Non-positive values indicate to use PETSc's default.
     */
    real_type minimum_step_size;

    /**
     * Maximum allowed step size for adaptive time stepping.
     * Non-positive values indicate to use PETSc's default.
     */
    real_type maximum_step_size;

    /**
     * Absolute error tolerance for adaptive time stepping.
     * Negative values indicate to use PETSc's default.
     */
    real_type absolute_tolerance;

    /**
     * Relative error tolerance for adaptive time stepping.
     * Negative values indicate to use PETSc's default.
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
   *  - a-la-PETSc style using TimeStepper::implicit_jacobian
   *    and TimeStepper::explicit_jacobian.
   *  - a-la-Deal.II style using TimeStepper::setup_jacobian and
   * TimeStepper::solve_for_jacobian_system
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
   * still using TimeStepper::solve_for_jacobian_system as a preconditioner.
   *
   * The first approach has instead the advantage that only the matrix assembly
   * procedure has to be coded, thus allowing quicker implementations and
   * faster turnaround for experimenting with linear solver preconditioning
   * configurations via command line customizations, like for example,
   * @code
   * ./myApp -ksp_type cg -pc_type gamg
   * @endcode
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
     * Change customization according to `data`.
     */
    void
    reinit(const TimeStepperData &data);

    /**
     * Set the preconditioning matrix only.
     *
     * When used with TimeStepper::setup_jacobian and
     * TimeStepper::solve_for_jacobian_system, PETSc will approximate the linear
     * system matrix-vector product using an internal matrix-free
     * representation.
     *
     * When used with TimeStepper::implicit_jacobian or
     * TimeStepper::explicit_jacobian, PETSc will use the same matrix for both
     * preconditioning and matrix-vector products.
     */
    void
    reinit_matrix(PMatrixType &P);

    /**
     * Set both the linear system matrix and the preconditioning matrix
     * that PETSc will use.
     */
    void
    reinit_matrices(AMatrixType &A, PMatrixType &P);

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
     * Here we also pass the matrix to handle the Jacobian.
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
     * Here we also pass the matrices to handle Jacobians.
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
     * y}$
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
     * Set up Jacobian callback without matrices.
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
     * Solvers must be provided via TimeStepper::solve_for_jacobian_system.
     */
    std::function<int(const real_type   t,
                      const VectorType &y,
                      const VectorType &ydot,
                      const real_type   alpha)>
      setup_jacobian;

    /**
     * Solution of the Jacobian system set up with TimeStepper::setup_jacobian.
     *
     * This is used as a preconditioner inside the Krylov process.
     */
    std::function<int(const VectorType &src, VectorType &dst)>
      solve_for_jacobian_system;

    /**
     * Return an index set containing the algebraic components.
     *
     * Implementation of this function is optional. If your equation is also
     * algebraic (i.e., it contains algebraic constraints, or Lagrange
     * multipliers), you should implement this function in order to return only
     * these components of your system.
     */
    std::function<IndexSet()> algebraic_components;

  protected:
    /**
     * The PETSc object
     */
    TS ts;

    SmartPointer<AMatrixType, TimeStepper> A;
    SmartPointer<PMatrixType, TimeStepper> P;
    bool                                   need_dae_tolerances;
  };

#  ifndef DOXYGEN
  /* Inline functions */

  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  inline TS
  TimeStepper<VectorType, PMatrixType, AMatrixType>::petsc_ts()
  {
    return ts;
  }

  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  inline MPI_Comm
  TimeStepper<VectorType, PMatrixType, AMatrixType>::get_mpi_communicator()
    const
  {
    return PetscObjectComm(reinterpret_cast<PetscObject>(ts));
  }

  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  inline typename TimeStepper<VectorType, PMatrixType, AMatrixType>::real_type
  TimeStepper<VectorType, PMatrixType, AMatrixType>::get_time()
  {
    PetscReal      t;
    PetscErrorCode ierr = TSGetTime(ts, &t);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    return t;
  }

  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  inline typename TimeStepper<VectorType, PMatrixType, AMatrixType>::real_type
  TimeStepper<VectorType, PMatrixType, AMatrixType>::get_time_step()
  {
    PetscReal      dt;
    PetscErrorCode ierr = TSGetTimeStep(ts, &dt);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    return dt;
  }
#  endif

} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

#endif
