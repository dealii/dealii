// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_sundials_arkode_h
#define dealii_sundials_arkode_h

#include <deal.II/base/config.h>


#ifdef DEAL_II_WITH_SUNDIALS

#  include <deal.II/base/conditional_ostream.h>
#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/mpi_stub.h>
#  include <deal.II/base/parameter_handler.h>

#  ifdef DEAL_II_WITH_PETSC
#    include <deal.II/lac/petsc_block_vector.h>
#    include <deal.II/lac/petsc_vector.h>
#  endif
#  include <deal.II/lac/vector.h>
#  include <deal.II/lac/vector_memory.h>

#  include <arkode/arkode.h>
#  include <nvector/nvector_serial.h>
#  ifdef DEAL_II_WITH_MPI
#    include <nvector/nvector_parallel.h>
#  endif
#  include <deal.II/base/discrete_time.h>

#  include <deal.II/sundials/arkode_stepper.h>
#  include <deal.II/sundials/n_vector.h>
#  include <deal.II/sundials/sundials_types.h>
#  include <deal.II/sundials/sunlinsol_wrapper.h>

#  include <boost/signals2.hpp>

#  include <sundials/sundials_linearsolver.h>
#  include <sundials/sundials_math.h>

#  include <exception>
#  include <memory>


DEAL_II_NAMESPACE_OPEN


/**
 * A namespace for dealing with ODE solvers through the SUNDIALS package.
 */
namespace SUNDIALS
{
  namespace internal
  {
    template <class Fn>
    struct FunctionProxy
    {
      std::function<Fn> *target;

      FunctionProxy(std::function<Fn> *t = nullptr)
        : target(t)
      {}

      FunctionProxy &
      operator=(std::function<Fn> *t)
      {
        target = t;
        return *this;
      }

      FunctionProxy &
      operator=(std::function<Fn> f)
      {
        Assert(
          target,
          ExcMessage(
            "Unable to assign the function since the target is not set in the "
            "proxy. Most probably, you tried to use it with a stepper that does "
            "not support this type of callback."));

        *target = std::move(f);
        return *this;
      }
    };
  } // namespace internal

  /**
   * Interface to SUNDIALS additive Runge-Kutta methods (ARKode).
   *
   * The class ARKode is a wrapper to SUNDIALS adaptive-step time integration
   * modules for stiff, nonstiff and mixed stiff/nonstiff systems of ordinary
   * differential equations (ODEs). This class provides access to general
   * underlying infrastructure of ARKODE, controls the main time marching loop,
   * output, reinitializaiton, exception handling. The class ARKode does not
   * implement any time integration strategies itself. Instead, it employs one
   * of the available stepper classes, which are basically wrappers around the
   * corresponding ARKODE steppers. Currently, only ARKStepper is available.
   * A user can also provide his own steppers or implement wrappers around,
   * e.g., ERKStepper, SPRKStep, LSRKStepper, or MRIStep modules of ARKode.
   *
   * In practice, you would of course at least want to set the end time up
   * to which the time integrator should run. In the example code shown here,
   * it is taken from the default value defined by ARKode::AdditionalData.
   */
  template <typename VectorType = Vector<double>>
  class ARKode
  {
  public:
    /**
     * Additional parameters that can be passed to the ARKode class.
     */
    class AdditionalData
    {
    public:
      /**
       * Initialization parameters for ARKode. Some parameters are provided for
       * the interface backward compatibility and have been effective moved to
       * the the AKRStepper additional parameters.
       *
       * Global parameters:
       *
       * @param initial_time Initial time
       * @param final_time Final time
       * @param initial_step_size Initial step size
       * @param output_period Desired time interval between each output
       *
       * Running parameters:
       *
       * @param minimum_step_size Minimum step size
       * @param maximum_order Maximum ARK order
       * @param maximum_non_linear_iterations Maximum number of nonlinear
       *   iterations
       * @param implicit_function_is_linear Specifies that the implicit portion
       *   of the problem is linear
       * @param implicit_function_is_time_independent Specifies that the
       *   implicit portion of the problem is linear and time independent
       * @param mass_is_time_independent Specifies that the mass pre-factor is
       *   independent of time
       * @param anderson_acceleration_subspace The number of vectors to use for
       *   Anderson acceleration within the packaged SUNDIALS solver.
       *
       * Error parameters:
       *
       * @param absolute_tolerance Absolute error tolerance
       * @param relative_tolerance Relative error tolerance
       */
      DEAL_II_DEPRECATED_WITH_COMMENT(
        "Use another constructor and ARKStepper::AdditionalData to initialize "
        "with the ARKStepper relevant parameters instead.")
      AdditionalData(
        // Initial parameters
        const double initial_time      = 0.0,
        const double final_time        = 1.0,
        const double initial_step_size = 1e-2,
        const double output_period     = 1e-1,
        // Running parameters
        const double       minimum_step_size                     = 1e-6,
        const unsigned int maximum_order                         = 5,
        const unsigned int maximum_non_linear_iterations         = 10,
        const bool         implicit_function_is_linear           = false,
        const bool         implicit_function_is_time_independent = false,
        const bool         mass_is_time_independent              = false,
        const int          anderson_acceleration_subspace        = 3,
        // Error parameters
        const double absolute_tolerance = 1e-6,
        const double relative_tolerance = 1e-5);

      /**
       * Initialization parameters for ARKode.
       *
       * Global parameters:
       *
       * @param initial_time Initial time
       * @param final_time Final time
       * @param initial_step_size Initial step size
       * @param output_period Desired time interval between each output
       *
       * Running parameters:
       *
       * @param minimum_step_size Minimum step size
       * @param maximum_order Maximum ARK order
       *
       * Error parameters:
       *
       * @param absolute_tolerance Absolute error tolerance
       * @param relative_tolerance Relative error tolerance
       */
      AdditionalData(
        // Initial parameters
        const double initial_time      = 0.0,
        const double final_time        = 1.0,
        const double initial_step_size = 1e-2,
        const double output_period     = 1e-1,
        // Running parameters
        const double       minimum_step_size = 1e-6,
        const unsigned int maximum_order     = 5,
        // Error parameters
        const double absolute_tolerance = 1e-6,
        const double relative_tolerance = 1e-5);

      /**
       * Add all AdditionalData() parameters to the given ParameterHandler
       * object. When the parameters are parsed from a file, the internal
       * parameters are automatically updated.
       *
       * The options you pass at construction time are set as default values in
       * the ParameterHandler object `prm`. You can later modify them by parsing
       * a parameter file using `prm`. The values of the parameter will be
       * updated whenever the content of `prm` is updated.
       *
       * Make sure that this class lives longer than `prm`. Undefined behavior
       * will occur if you destroy this class, and then parse a parameter file
       * using `prm`.
       */
      void
      add_parameters(ParameterHandler &prm);

      /**
       * Initial time for the DAE.
       */
      double initial_time;

      /**
       * Final time.
       */
      double final_time;

      /**
       * Initial step size.
       */
      double initial_step_size;

      /**
       * Minimum step size.
       */
      double minimum_step_size;

      /**
       * Absolute error tolerance for adaptive time stepping.
       */
      double absolute_tolerance;

      /**
       * Relative error tolerance for adaptive time stepping.
       */
      double relative_tolerance;

      /**
       * Maximum order of ARK.
       */
      unsigned int maximum_order;

      /**
       * Desired time period between each output. The actual output time period
       * may be adjusted by Arkode.
       */
      double output_period;

      /**
       * Maximum number of iterations for Newton or fixed point method during
       * time advancement.
       *
       * @deprecated Specify ARKStepper::AdditionalData::maximum_non_linear_iterations instead.
       */
      DEAL_II_DEPRECATED_WITH_COMMENT(
        "Specify ARKStepper::AdditionalData::maximum_non_linear_iterations "
        "instead.")
      unsigned int maximum_non_linear_iterations;

      /**
       * Specify whether the implicit portion of the problem is linear.
       *
       * @deprecated Specify ARKStepper::AdditionalData::implicit_function_is_linear instead.
       */
      DEAL_II_DEPRECATED_WITH_COMMENT(
        "Specify ARKStepper::AdditionalData::implicit_function_is_linear "
        "instead.")
      bool implicit_function_is_linear;

      /**
       * Specify whether the implicit portion of the problem is linear and time
       * independent.
       *
       * @deprecated Specify ARKStepper::AdditionalData::implicit_function_is_time_independent instead.
       */
      DEAL_II_DEPRECATED_WITH_COMMENT(
        "Specify ARKStepper::AdditionalData::implicit_function_is_time_independent "
        "instead.")
      bool implicit_function_is_time_independent;

      /**
       * Specify whether the mass pre-factor is time independent. Has no effect
       * if no mass is specified.
       *
       * @deprecated Specify ARKStepper::AdditionalData::mass_is_time_independent instead.
       */
      DEAL_II_DEPRECATED_WITH_COMMENT(
        "Specify ARKStepper::AdditionalData::mass_is_time_independent "
        "instead.")
      bool mass_is_time_independent;

      /**
       * Number of subspace vectors to use for Anderson acceleration. Only
       * meaningful if the packaged SUNDIALS fixed-point solver is used.
       *
       * @deprecated Specify ARKStepper::AdditionalData::anderson_acceleration_subspace instead.
       */
      DEAL_II_DEPRECATED_WITH_COMMENT(
        "Specify ARKStepper::AdditionalData::anderson_acceleration_subspace "
        "instead.")
      int anderson_acceleration_subspace;
    };

    /**
     * Constructor, with class parameters set by the AdditionalData object.
     *
     * @param data ARKode configuration data
     *
     * @note With SUNDIALS 6 and later this constructor sets up logging
     * objects to only work on the present processor (i.e., results are only
     * communicated over MPI_COMM_SELF).
     */
    ARKode(const AdditionalData &data = AdditionalData());

    /**
     * Constructor receiving the additional data and MPI communicator.
     *
     * @param data ARKode configuration data
     * @param mpi_comm MPI Communicator over which logging operations are
     * computed. Only used in SUNDIALS 6 and newer.
     */
    ARKode(const AdditionalData &data, const MPI_Comm mpi_comm);

    /**
     * Constructor receiving the stepper and optionally the additional data.
     *
     * @param stepper ARKodeStepper object
     * @param data ARKode configuration data
     *
     * @note With SUNDIALS 6 and later this constructor sets up logging
     * objects to only work on the present processor (i.e., results are only
     * communicated over MPI_COMM_SELF).
     */
    ARKode(std::shared_ptr<ARKodeStepper<VectorType>> stepper,
           const AdditionalData                      &data = AdditionalData());

    /**
     * Constructor receiving the stepper, additional data and MPI communicator.
     *
     * @param stepper ARKodeStepper object
     * @param data ARKode configuration data
     * @param mpi_comm MPI Communicator over which logging operations are
     * computed. Only used in SUNDIALS 6 and newer.
     */
    ARKode(std::shared_ptr<ARKodeStepper<VectorType>> stepper,
           const AdditionalData                      &data,
           const MPI_Comm                             mpi_comm);

    /**
     * Destructor.
     */
    ~ARKode();

    /**
     * Integrate the initial value problem. This function returns the final
     * number of computed steps.
     *
     * @param solution On input, this vector contains the initial condition. On
     *   output, it contains the solution at the final time.
     */
    unsigned int
    solve_ode(VectorType &solution);

    /**
     * Integrate the initial value problem. Compared to the function above, this
     * function allows to specify an @p intermediate_time for the next solution.
     * Repeated calls of this function must use monotonously increasing values
     * for @p intermediate_time. The last solution state is saved internally
     * along with the @p intermediate_time and will be reused as initial
     * condition for the next call.
     *
     * Users may find this function useful when integrating ARKode into an outer
     * time loop of their own, especially when output_step() is too restrictive.
     *
     * @note @p intermediate_time may be larger than AdditionalData::final_time,
     *   which is ignored by this function.
     *
     * @param solution The final solution. If the solver restarts, either
     *   because it is the first ever solve or the flag @p reset_solver is
     *   set, the vector is also used as initial condition.
     * @param intermediate_time The time for the incremental solution step. Must
     *   be greater than the last time that was used in a previous call to this
     *   function.
     * @param reset_solver Optional flag to recreate all internal objects which
     *   may be desirable for spatial adaptivity methods. If set to `true`,
     *   reset() is called before solving the ODE, which sets @p solution as
     *   initial condition. This will *not* reset the stored time from previous
     *   calls to this function.
     */
    unsigned int
    solve_ode_incrementally(VectorType  &solution,
                            const double intermediate_time,
                            const bool   reset_solver = false);

    /**
     * Clear internal memory and start with clean objects. This function is
     * called when the simulation starts and when the user returns true to a
     * call to solver_should_restart().
     *
     * By default solver_should_restart() returns false. If the user needs to
     * implement, for example, local adaptivity in space, he or she may assign
     * a different function to solver_should_restart() that performs all mesh
     * changes, transfers the solution to the new mesh, and returns true.
     *
     * @param t  The new starting time
     * @param h  The new starting time step
     * @param y  The new initial solution
     */
    void
    reset(const double t, const double h, const VectorType &y);

    /**
     * Provides user access to the internally used ARKODE memory.
     *
     * This functionality is intended for users who wish to query additional
     * information directly from the ARKODE integrator, refer to the ARKODE
     * manual for the various `ARKStepGet...` functions. The `ARKStepSet...`
     * functions should not be called since this might lead to conflicts with
     * various settings that are performed by this ARKode object.
     *
     * @note If custom settings of ARKODE functionality (that are not achievable
     *   via the interface of this class) are required, the function
     *   custom_setup() should be used.
     *
     * @return pointer to the ARKODE memory block that can be passed to SUNDIALS
     *   functions
     */
    void *
    get_arkode_memory() const;

    /**
     * A function object that users may supply and that is intended to compute
     * the explicit part of the IVP right hand side. Sets $explicit_f = f_E(t,
     * y)$.
     *
     * @note This member represents a local proxy wrapper around the corresponding
     * property of ARKStepper and should be used only if ARKStepper is provided
     * to ARKode as a stepper. If this requirement is violated, an exception is
     * triggered upon access to this field.
     *
     * @deprecated Specify @p ARKStepper::explicit_function instead.
     */
    internal::FunctionProxy<
      void(const double t, const VectorType &y, VectorType &explicit_f)>
      explicit_function;

    /**
     * A function object that users may supply and that is intended to compute
     * the implicit part of the IVP right hand side. Sets $implicit_f = f_I(t,
     * y)$.
     *
     * @note This member represents a local proxy wrapper around the corresponding
     * property of ARKStepper and should be used only if ARKStepper is provided
     * to ARKode as a stepper. If this requirement is violated, an exception is
     * triggered upon access to this field.
     *
     * @deprecated Specify @p ARKStepper::implicit_function instead.
     */
    internal::FunctionProxy<
      void(const double t, const VectorType &y, VectorType &res)>
      implicit_function;

    /**
     * A function object that users may supply and that is intended to compute
     * the product of the @ref GlossMassMatrix "mass matrix" with a given vector `v`. This function
     * will be called by ARKode (possibly several times) after
     * mass_times_setup() has been called at least once. ARKode tries to do its
     * best to call mass_times_setup() the minimum amount of times.
     *
     * @note This member represents a local proxy wrapper around the corresponding
     * property of ARKStepper and should be used only if ARKStepper is provided
     * to ARKode as a stepper. If this requirement is violated, an exception is
     * triggered upon access to this field.
     *
     * @deprecated Specify @p ARKStepper::mass_times_vector instead.
     */
    internal::FunctionProxy<
      void(const double t, const VectorType &v, VectorType &Mv)>
      mass_times_vector;

    /**
     * A function object that users may supply and that is intended to set up
     * the @ref GlossMassMatrix "mass matrix". This function is called by ARKode any time a mass
     * matrix update is required. The user should compute the mass matrix (or
     * update all the variables that allow the application of the mass matrix).
     * This function is guaranteed to be called by ARKode at least once, before
     * any call to mass_times_vector().
     *
     * @note This member represents a local proxy wrapper around the corresponding
     * property of ARKStepper and should be used only if ARKStepper is provided
     * to ARKode as a stepper. If this requirement is violated, an exception is
     * triggered upon access to this field.
     *
     * @deprecated Specify @p ARKStepper::mass_times_setup instead.
     */
    internal::FunctionProxy<void(const double t)> mass_times_setup;

    /**
     * A function object that users may supply and that is intended to compute
     * the product of the Jacobian matrix with a given vector `v`. The Jacobian
     * here refers to $J=\frac{\partial f_I}{\partial y}$, i.e., the Jacobian of
     * the user-specified implicit_function.
     *
     * @note This member represents a local proxy wrapper around the corresponding
     * property of ARKStepper and should be used only if ARKStepper is provided
     * to ARKode as a stepper. If this requirement is violated, an exception is
     * triggered upon access to this field.
     *
     * @deprecated Specify @p ARKStepper::jacobian_times_vector instead.
     */
    internal::FunctionProxy<void(const VectorType &v,
                                 VectorType       &Jv,
                                 const double      t,
                                 const VectorType &y,
                                 const VectorType &fy)>
      jacobian_times_vector;

    /**
     * A function object that users may supply and that is intended to set up
     * all data necessary for the application of jacobian_times_vector().
     *
     * @note This member represents a local proxy wrapper around the corresponding
     * property of ARKStepper and should be used only if ARKStepper is provided
     * to ARKode as a stepper. If this requirement is violated, an exception is
     * triggered upon access to this field.
     *
     * @deprecated Specify @p ARKStepper::jacobian_times_setup instead.
     */
    internal::FunctionProxy<
      void(const double t, const VectorType &y, const VectorType &fy)>
      jacobian_times_setup;

    /**
     * A LinearSolveFunction object that users may supply and that is intended
     * to solve the linearized system $Ax=b$, where $A = M-\gamma J$ is the
     * Jacobian of the nonlinear residual.
     *
     * @note This member represents a local proxy wrapper around the corresponding
     * property of ARKStepper and should be used only if ARKStepper is provided
     * to ARKode as a stepper. If this requirement is violated, an exception is
     * triggered upon access to this field.
     *
     * @deprecated Specify @p ARKStepper::solve_linearized_system instead.
     */
    internal::FunctionProxy<void(SundialsOperator<VectorType>       &op,
                                 SundialsPreconditioner<VectorType> &prec,
                                 VectorType                         &x,
                                 const VectorType                   &b,
                                 double                              tol)>
      solve_linearized_system;

    /**
     * A LinearSolveFunction object that users may supply and that is intended
     * to solve the mass system $Mx=b$.
     *
     * @note This member represents a local proxy wrapper around the corresponding
     * property of ARKStepper and should be used only if ARKStepper is provided
     * to ARKode as a stepper. If this requirement is violated, an exception is
     * triggered upon access to this field.
     *
     * @deprecated Specify @p ARKStepper::solve_mass instead.
     */
    internal::FunctionProxy<void(SundialsOperator<VectorType>       &op,
                                 SundialsPreconditioner<VectorType> &prec,
                                 VectorType                         &x,
                                 const VectorType                   &b,
                                 double                              tol)>
      solve_mass;

    /**
     * A function object that users may supply to either pass a preconditioner
     * to a SUNDIALS built-in solver or to apply a custom preconditioner within
     * the user's own linear solve specified in solve_linearized_system().
     *
     * @note This member represents a local proxy wrapper around the corresponding
     * property of ARKStepper and should be used only if ARKStepper is provided
     * to ARKode as a stepper. If this requirement is violated, an exception is
     * triggered upon access to this field.
     *
     * @deprecated Specify @p ARKStepper::jacobian_preconditioner_solve instead.
     */
    internal::FunctionProxy<void(const double      t,
                                 const VectorType &y,
                                 const VectorType &fy,
                                 const VectorType &r,
                                 VectorType       &z,
                                 const double      gamma,
                                 const double      tol,
                                 const int         lr)>
      jacobian_preconditioner_solve;

    /**
     * A function object that users may supply to set up a preconditioner
     * specified in jacobian_preconditioner_solve().
     *
     * @note This member represents a local proxy wrapper around the corresponding
     * property of ARKStepper and should be used only if ARKStepper is provided
     * to ARKode as a stepper. If this requirement is violated, an exception is
     * triggered upon access to this field.
     *
     * @deprecated Specify @p ARKStepper::jacobian_preconditioner_setup instead.
     */
    internal::FunctionProxy<void(const double      t,
                                 const VectorType &y,
                                 const VectorType &fy,
                                 const int         jok,
                                 int              &jcur,
                                 const double      gamma)>
      jacobian_preconditioner_setup;

    /**
     * A function object that users may supply to either pass a preconditioner
     * to a SUNDIALS built-in solver or to apply a custom preconditioner within
     * the user's own linear solve specified in solve_mass().
     *
     * @note This member represents a local proxy wrapper around the corresponding
     * property of ARKStepper and should be used only if ARKStepper is provided
     * to ARKode as a stepper. If this requirement is violated, an exception is
     * triggered upon access to this field.
     *
     * @deprecated Specify @p ARKStepper::mass_preconditioner_solve instead.
     */
    internal::FunctionProxy<void(const double      t,
                                 const VectorType &r,
                                 VectorType       &z,
                                 const double      tol,
                                 const int         lr)>
      mass_preconditioner_solve;

    /**
     * A function object that users may supply to set up a preconditioner
     * specified in mass_preconditioner_setup().
     *
     * @note This member represents a local proxy wrapper around the corresponding
     * property of ARKStepper and should be used only if ARKStepper is provided
     * to ARKode as a stepper. If this requirement is violated, an exception is
     * triggered upon access to this field.
     *
     * @deprecated Specify @p ARKStepper::mass_preconditioner_solve instead.
     */
    internal::FunctionProxy<void(const double t)> mass_preconditioner_setup;

    /**
     * A function object that users may supply and that is intended to
     * postprocess the solution. This function is called by ARKode at fixed
     * time increments (every `output_period` time units), and it is passed a
     * polynomial interpolation of the solution, computed using the current ARK
     * order and the (internally stored) previously computed solution steps.
     *
     * @note It is well possible that internally ARKode computes a time
     *   step which is much larger than the `output_period` step, and therefore
     *   calls this function consecutively several times by simply performing
     *   all intermediate interpolations. There is no relationship between how
     *   many times this function is called and how many time steps have
     *   actually been computed.
     */
    std::function<void(const double       t,
                       const VectorType  &sol,
                       const unsigned int step_number)>
      output_step;

    /**
     * A function object that users may supply and that is intended to evaluate
     * whether the solver should be restarted (for example because the number of
     * degrees of freedom has changed).
     *
     * This function is supposed to perform all operations that are necessary
     * in `sol` to make sure that the resulting vectors are consistent, and of
     * the correct final size.
     *
     * For example, one may decide that a local refinement is necessary at time
     * t. This function should then return true, and change the dimension of
     * `sol` to reflect the new dimension. Since ARKode does not know about the
     * new dimension, an internal reset is necessary.
     *
     * The default implementation simply returns `false`, i.e., no restart is
     * performed during the evolution.
     */
    std::function<bool(const double t, VectorType &sol)> solver_should_restart;

    /**
     * A function object that users may supply and that is intended to return a
     * vector whose components are the weights used by ARKode to compute the
     * vector norm. The implementation of this function is optional, and it is
     * used only if implemented.
     */
    std::function<VectorType &()> get_local_tolerances;

    /**
     * A function object that users may supply and which is intended to perform
     * custom settings on the supplied @p arkode_mem object. Refer to the
     * SUNDIALS documentation for valid options.
     *
     * For instance, the following code sets Lagrange interpolants to be used
     * for dense output (interpolation of solution output values) and implicit
     * method predictors:
     *
     * @code
     *      ode.custom_setup = [&](void *arkode_mem) {
     *        ARKodeSetInterpolantType(arkode_mem, ARK_INTERP_LAGRANGE);
     *      };
     * @endcode
     *
     * @note This function will be called at the end of all other set up right
     *   before the actual time evolution is started or continued with
     *   solve_ode(). This function is also called when the solver is restarted,
     *   see solver_should_restart(). Consult the SUNDIALS manual to see which
     *   options are still available at this point.
     *
     * @param arkode_mem pointer to the ARKODE memory block which can be used
     *   for custom calls to `ARKodeSet...` methods.
     */
    std::function<void(void *arkode_mem)> custom_setup;

  private:
    /**
     * Throw an exception when a function with the given name is not
     * implemented.
     */
    DeclException1(ExcFunctionNotProvided,
                   std::string,
                   << "Please provide an implementation for the function \""
                   << arg1 << "\"");

    /**
     * Internal routine to call ARKode repeatedly.
     */
    unsigned int
    do_evolve_time(VectorType           &solution,
                   dealii::DiscreteTime &time,
                   const bool            do_reset);

    /**
     * This function is executed at construction time to set the
     * std::function above to trigger an assert if they are not
     * implemented.
     */
    void
    set_functions_to_trigger_an_assert();

    /**
     * ARKode configuration data.
     */
    AdditionalData data;

    /**
     * Stepper object.
     */
    std::shared_ptr<ARKodeStepper<VectorType>> stepper;

#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    /**
     * A context object associated with the ARKode solver.
     */
    SUNContext arkode_ctx;
#  endif

    /**
     * MPI communicator. Only used for SUNDIALS' logging routines - the actual
     * solve routines will use the communicator provided by the vector class.
     */
    MPI_Comm mpi_communicator;

    /**
     * The final time in the last call to solve_ode().
     */
    double last_end_time;

    /**
     * A pointer to any exception that may have been thrown in user-defined
     * call-backs and that we have to deal after the KINSOL function we call
     * has returned.
     */
    mutable std::exception_ptr pending_exception;

#  ifdef DEAL_II_WITH_PETSC
#    ifdef PETSC_USE_COMPLEX
    static_assert(!std::is_same_v<VectorType, PETScWrappers::MPI::Vector>,
                  "Sundials does not support complex scalar types, "
                  "but PETSc is configured to use a complex scalar type!");

    static_assert(!std::is_same_v<VectorType, PETScWrappers::MPI::BlockVector>,
                  "Sundials does not support complex scalar types, "
                  "but PETSc is configured to use a complex scalar type!");
#    endif // PETSC_USE_COMPLEX
#  endif   // DEAL_II_WITH_PETSC
  };


  template <typename VectorType>
  ARKode<VectorType>::AdditionalData::AdditionalData(
    // Initial parameters
    const double initial_time,
    const double final_time,
    const double initial_step_size,
    const double output_period,
    // Running parameters
    const double       minimum_step_size,
    const unsigned int maximum_order,
    const unsigned int maximum_non_linear_iterations,
    const bool         implicit_function_is_linear,
    const bool         implicit_function_is_time_independent,
    const bool         mass_is_time_independent,
    const int          anderson_acceleration_subspace,
    // Error parameters
    const double absolute_tolerance,
    const double relative_tolerance)
    : initial_time(initial_time)
    , final_time(final_time)
    , initial_step_size(initial_step_size)
    , minimum_step_size(minimum_step_size)
    , absolute_tolerance(absolute_tolerance)
    , relative_tolerance(relative_tolerance)
    , maximum_order(maximum_order)
    , output_period(output_period)
    , maximum_non_linear_iterations(maximum_non_linear_iterations)
    , implicit_function_is_linear(implicit_function_is_linear)
    , implicit_function_is_time_independent(
        implicit_function_is_time_independent)
    , mass_is_time_independent(mass_is_time_independent)
    , anderson_acceleration_subspace(anderson_acceleration_subspace)
  {}

  template <typename VectorType>
  ARKode<VectorType>::AdditionalData::AdditionalData(
    // Initial parameters
    const double initial_time,
    const double final_time,
    const double initial_step_size,
    const double output_period,
    // Running parameters
    const double       minimum_step_size,
    const unsigned int maximum_order,
    // Error parameters
    const double absolute_tolerance,
    const double relative_tolerance)
    : initial_time(initial_time)
    , final_time(final_time)
    , initial_step_size(initial_step_size)
    , minimum_step_size(minimum_step_size)
    , absolute_tolerance(absolute_tolerance)
    , relative_tolerance(relative_tolerance)
    , maximum_order(maximum_order)
    , output_period(output_period)
    , maximum_non_linear_iterations(10)
    , implicit_function_is_linear(false)
    , implicit_function_is_time_independent(false)
    , mass_is_time_independent(false)
    , anderson_acceleration_subspace(3)
  {}

  template <typename VectorType>
  void
  ARKode<VectorType>::AdditionalData::add_parameters(ParameterHandler &prm)
  {
    prm.add_parameter("Initial time", initial_time);
    prm.add_parameter("Final time", final_time);
    prm.add_parameter("Time interval between each output", output_period);
    prm.enter_subsection("Running parameters");
    prm.add_parameter("Initial step size", initial_step_size);
    prm.add_parameter("Minimum step size", minimum_step_size);
    prm.add_parameter("Maximum order of ARK", maximum_order);
    prm.add_parameter("Maximum number of nonlinear iterations",
                      maximum_non_linear_iterations);
    prm.add_parameter("Implicit function is linear",
                      implicit_function_is_linear);
    prm.add_parameter("Implicit function is time independent",
                      implicit_function_is_time_independent);
    prm.add_parameter("Mass is time independent", mass_is_time_independent);
    prm.add_parameter("Anderson-acceleration subspace",
                      anderson_acceleration_subspace);
    prm.leave_subsection();
    prm.enter_subsection("Error control");
    prm.add_parameter("Absolute error tolerance", absolute_tolerance);
    prm.add_parameter("Relative error tolerance", relative_tolerance);
    prm.leave_subsection();
  }

} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif


#endif
