// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_sundials_arkode_stepper_h
#define dealii_sundials_arkode_stepper_h

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

#  include <deal.II/sundials/invocation_context.h>
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
    template <typename Stepper>
    struct ARKCallbackContext
    {
      Stepper            *stepper           = nullptr;
      std::exception_ptr *pending_exception = nullptr;
    };
  } // namespace internal

  template <typename VectorType>
  class ARKodeStepper
  {
  public:
    virtual ~ARKodeStepper() = default;

    virtual void
    reinit(double                      t0,
           const VectorType           &y0,
           internal::InvocationContext inv_ctx) = 0;

    /**
     * Provides user access to the internally used ARKODE memory.
     *
     * This functionality is intended for users who wish to query additional
     * information directly from the ARKODE integrator, refer to the ARKODE
     * manual for the various `ARKodeGet...` functions. The `ARKodeSet...`
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
    virtual void *
    get_arkode_memory() const = 0;
  };


  template <typename VectorType = Vector<double>>
  class ARKStepper : public ARKodeStepper<VectorType>
  {
  public:
    /**
     * Additional parameters that can be passed to the ARKStepper class.
     */
    class AdditionalData
    {
    public:
      /**
       * Initialization parameters for ARKStepper.
       *
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
       */
      AdditionalData(const unsigned int maximum_non_linear_iterations = 10,
                     const bool         implicit_function_is_linear   = false,
                     const bool implicit_function_is_time_independent = false,
                     const bool mass_is_time_independent              = false,
                     const int  anderson_acceleration_subspace        = 3);

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
       * Maximum number of iterations for Newton or fixed point method during
       * time advancement.
       */
      unsigned int maximum_non_linear_iterations;

      /**
       * Specify whether the implicit portion of the problem is linear.
       */
      bool implicit_function_is_linear;

      /**
       * Specify whether the implicit portion of the problem is linear and time
       * independent.
       */
      bool implicit_function_is_time_independent;

      /**
       * Specify whether the mass pre-factor is time independent. Has no effect
       * if no mass is specified.
       */
      bool mass_is_time_independent;

      /**
       * Number of subspace vectors to use for Anderson acceleration. Only
       * meaningful if the packaged SUNDIALS fixed-point solver is used.
       */
      int anderson_acceleration_subspace;
    };

    ARKStepper(const AdditionalData &data = AdditionalData());

    ~ARKStepper();

    void
    reinit(double                      t0,
           const VectorType           &y0,
           internal::InvocationContext inv_ctx) override;

    void *
    get_arkode_memory() const override;

    /**
     * A function object that users may supply and that is intended to compute
     * the explicit part of the IVP right hand side. Sets $explicit_f = f_E(t,
     * y)$.
     *
     * At least one of explicit_function() or implicit_function() must be
     * provided. According to which one is provided, explicit, implicit, or
     * mixed RK methods are used.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, ARKode can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
     */
    std::function<
      void(const double t, const VectorType &y, VectorType &explicit_f)>
      explicit_function;

    /**
     * A function object that users may supply and that is intended to compute
     * the implicit part of the IVP right hand side. Sets $implicit_f = f_I(t,
     * y)$.
     *
     * At least one of explicit_function() or implicit_function() must be
     * provided. According to which one is provided, explicit, implicit, or
     * mixed RK methods are used.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, ARKode can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
     */
    std::function<void(const double t, const VectorType &y, VectorType &res)>
      implicit_function;

    /**
     * A function object that users may supply and that is intended to compute
     * the product of the @ref GlossMassMatrix "mass matrix" with a given vector `v`. This function
     * will be called by ARKode (possibly several times) after
     * mass_times_setup() has been called at least once. ARKode tries to do its
     * best to call mass_times_setup() the minimum amount of times.
     *
     * A call to this function should store in `Mv` the result of $M$
     * applied to `v`.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, ARKode can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
     */
    std::function<void(const double t, const VectorType &v, VectorType &Mv)>
      mass_times_vector;

    /**
     * A function object that users may supply and that is intended to set up
     * the @ref GlossMassMatrix "mass matrix". This function is called by ARKode any time a mass
     * matrix update is required. The user should compute the mass matrix (or
     * update all the variables that allow the application of the mass matrix).
     * This function is guaranteed to be called by ARKode at least once, before
     * any call to mass_times_vector().
     *
     * ARKode supports the case where the mass matrix may depend on time, but
     * not the case where the mass matrix depends on the solution itself.
     *
     * If the user does not provide a mass_times_vector() function, then the
     * identity is used. If the mass_times_setup() function is not provided,
     * then mass_times_vector() should do all the work by itself.
     *
     * If the user uses a matrix-based computation of the mass matrix, then
     * this is the right place where an assembly routine should be called to
     * assemble the matrix. Subsequent calls (possibly  more than one) to
     * mass_times_vector() can assume that this function has been called at
     * least once.
     *
     * @note No assumption is made by this interface on what the user
     *   should do in this function. ARKode only assumes that after a call to
     *   mass_times_setup() it is possible to call mass_times_vector().
     *
     * @param t The current evaluation time
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, ARKode can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
     */
    std::function<void(const double t)> mass_times_setup;

    /**
     * A function object that users may supply and that is intended to compute
     * the product of the Jacobian matrix with a given vector `v`. The Jacobian
     * here refers to $J=\frac{\partial f_I}{\partial y}$, i.e., the Jacobian of
     * the user-specified implicit_function.
     *
     * A call to this function should store in `Jv` the result of $J$
     * applied to `v`.
     *
     * Arguments to the function are
     *
     * @param[in] v  The vector to be multiplied by the Jacobian
     * @param[out] Jv The vector to be filled with the product J*v
     * @param[in] t  The current time
     * @param[in] y  The current $y$ vector for the current ARKode internal
     *   step
     * @param[in] fy  The current value of the implicit right-hand side at y,
     *   $f_I (t_n, y)$.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, ARKode can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
     */
    std::function<void(const VectorType &v,
                       VectorType       &Jv,
                       const double      t,
                       const VectorType &y,
                       const VectorType &fy)>
      jacobian_times_vector;

    /**
     * A function object that users may supply and that is intended to set up
     * all data necessary for the application of jacobian_times_vector(). This
     * function is called by ARKode any time a Jacobian update is required.
     * The user should compute the Jacobian (or update all the variables that
     * allow the application of Jacobian). This function is guaranteed to
     * be called by ARKode at least once, before any call to
     * jacobian_times_vector().
     *
     * If the jacobian_times_setup() function is not provided, then
     * jacobian_times_vector() should do all the work by itself.
     *
     * If the user uses a matrix based computation of the Jacobian, then this is
     * the right place where an assembly routine should be called to assemble
     * the matrix. Subsequent calls (possibly  more than one) to
     * jacobian_times_vector() can assume that this function has been called at
     * least once.
     *
     * @note No assumption is made by this interface on what the user
     *   should do in this function. ARKode only assumes that after a call to
     *   jacobian_times_setup() it is possible to call jacobian_times_vector().
     *
     * @param t  The current time
     * @param y  The current ARKode internal solution vector $y$
     * @param fy  The implicit right-hand side function evaluated at the
     *   current time $t$ and state $y$, i.e., $f_I(y,t)$
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, ARKode can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
     */
    std::function<
      void(const double t, const VectorType &y, const VectorType &fy)>
      jacobian_times_setup;

    /**
     * A LinearSolveFunction object that users may supply and that is intended
     * to solve the linearized system $Ax=b$, where $A = M-\gamma J$ is the
     * Jacobian of the nonlinear residual. The application of the @ref GlossMassMatrix "mass matrix"
     * $M$ and Jacobian $J$ are known through the functions mass_times_vector()
     * and jacobian_times_vector() and $\gamma$ is a factor provided by
     * SUNDIALS. The matrix-vector product $Ax$ is encoded in the supplied
     * SundialsOperator. If a preconditioner was set through
     * jacobian_preconditioner_solve(), it is encoded in the
     * SundialsPreconditioner. If no preconditioner was supplied this way, the
     * preconditioner is the identity matrix, i.e., no preconditioner. The user
     * is free to use a custom preconditioner in this function object that is
     * not supplied through SUNDIALS.
     *
     * If you do not specify a solve_linearized_system() function, then a
     * SUNDIALS packaged SPGMR solver with default settings is used.
     *
     * For more details on the function type refer to LinearSolveFunction.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, ARKode can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
     */
    LinearSolveFunction<VectorType> solve_linearized_system;

    /**
     * A LinearSolveFunction object that users may supply and that is intended
     * to solve the mass system $Mx=b$. The matrix-vector product $Mx$ is
     * encoded in the supplied SundialsOperator. If a preconditioner was set
     * through mass_preconditioner_solve(), it is encoded in the
     * SundialsPreconditioner. If no preconditioner was supplied this way, the
     * preconditioner is the identity matrix, i.e., no preconditioner. The user
     * is free to use a custom preconditioner in this function object that is
     * not supplied through SUNDIALS.
     *
     * The user must specify this function if a non-identity mass matrix is used
     * and applied in mass_times_vector().
     *
     * For more details on the function type refer to LinearSolveFunction.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, ARKode can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
     */
    LinearSolveFunction<VectorType> solve_mass;

    /**
     * A function object that users may supply to either pass a preconditioner
     * to a SUNDIALS built-in solver or to apply a custom preconditioner within
     * the user's own linear solve specified in solve_linearized_system().
     *
     * This function should compute the solution to the preconditioner equation
     * $Pz=r$ and store it in @p z. In this equation $P$ should approximate the
     * Jacobian $M-\gamma J$ of the nonlinear system.
     *
     * @param[in] t  The current time
     * @param[in] y  The current $y$ vector for the current ARKode internal
     *   step
     * @param[in] fy  The current value of the implicit right-hand side at y,
     *   $f_I (t_n, y)$.
     * @param[in] r  The right-hand side of the preconditioner equation
     * @param[out] z The solution of applying the preconditioner, i.e., solving
     *   $Pz=r$
     * @param[in] gamma The value $\gamma$ in the preconditioner equation
     * @param[in] tol The tolerance up to which the system should be solved
     * @param[in] lr An input flag indicating whether the preconditioner solve
     *   is to use the left preconditioner (lr = 1) or the right preconditioner
     *   (lr = 2). Only relevant if used with a SUNDIALS packaged solver. If
     *   used with a custom solve_mass() function this will be set to zero.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, ARKode can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
     */
    std::function<void(const double      t,
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
     * This function should prepare the solution of the preconditioner equation
     * $Pz=r$. In this equation $P$ should approximate the Jacobian $M-\gamma J$
     * of the nonlinear system.
     *
     * If the jacobian_preconditioner_setup() function is not provided, then
     * jacobian_preconditioner_solve() should do all the work by itself.
     *
     * @note No assumption is made by this interface on what the user
     *   should do in this function. ARKode only assumes that after a call to
     *   jacobian_preconditioner_setup() it is possible to call
     *   jacobian_preconditioner_solve().
     *
     * @param[in] t  The current time
     * @param[in] y  The current $y$ vector for the current ARKode internal
     *   step
     * @param[in] fy  The current value of the implicit right-hand side at y,
     *   $f_I (t_n, y)$.
     * @param[in] jok  An input flag indicating whether the Jacobian-related
     *   data needs to be updated. The jok argument provides for the reuse of
     *   Jacobian data in the preconditioner solve function. When jok =
     *   SUNFALSE, the Jacobian-related data should be recomputed from scratch.
     *   When jok = SUNTRUE the Jacobian data, if saved from the previous call
     *   to this function, can be reused (with the current value of gamma). A
     *   call with jok = SUNTRUE can only occur after a call with jok =
     *   SUNFALSE.
     * @param[out] jcur On output this should be set to SUNTRUE if Jacobian data
     *   was recomputed, or set to SUNFALSE if Jacobian data was not recomputed,
     *   but saved data was still reused.
     * @param[in] gamma The value $\gamma$ in $M-\gamma J$. The preconditioner
     *   should approximate the inverse of this matrix.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, ARKode can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
     */
    std::function<void(const double      t,
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
     * This function should compute the solution to the preconditioner equation
     * $Pz=r$ and store it in @p z. In this equation $P$ should approximate the
     * @ref GlossMassMatrix "mass matrix" $M$.
     *
     * @param[in] t  The current time
     * @param[in] r  The right-hand side of the preconditioner equation
     * @param[out] z The solution of applying the preconditioner, i.e., solving
     *   $Pz=r$
     * @param[in] gamma The value $\gamma$ in the preconditioner equation
     * @param[in] tol The tolerance up to which the system should be solved
     * @param[in] lr An input flag indicating whether the preconditioner solve
     *   is to use the left preconditioner (lr = 1) or the right preconditioner
     *   (lr = 2). Only relevant if used with a SUNDIALS packaged solver. If
     *   used with a custom solve_mass() function this will be set to zero.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, ARKode can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
     */
    std::function<void(const double      t,
                       const VectorType &r,
                       VectorType       &z,
                       const double      tol,
                       const int         lr)>
      mass_preconditioner_solve;

    /**
     * A function object that users may supply to set up a preconditioner
     * specified in mass_preconditioner_solve().
     *
     * This function should prepare the solution of the preconditioner equation
     * $Pz=r$. In this equation $P$ should approximate the @ref GlossMassMatrix "mass matrix" $M$.
     *
     * If the mass_preconditioner_setup() function is not provided, then
     * mass_preconditioner_solve() should do all the work by itself.
     *
     * @note No assumption is made by this interface on what the user
     *   should do in this function. ARKode only assumes that after a call to
     *   mass_preconditioner_setup() it is possible to call
     *   mass_preconditioner_solve().
     *
     * @param[in] t  The current time
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, ARKode can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
     */
    std::function<void(const double t)> mass_preconditioner_setup;

    /**
     * A function object that users may supply and which is intended to perform
     * custom settings on the supplied @p arkode_mem object. Refer to the
     * SUNDIALS documentation for valid options.
     *
     * For instance, the following code attaches two files for diagnostic and
     * error output of the internal ARKODE implementation:
     *
     * @code
     *      ode.custom_setup = [&](void *arkode_mem) {
     *        ARKStepSetErrFile(arkode_mem, errfile);
     *        ARKStepSetDiagnostics(arkode_mem, diagnostics_file);
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
     *   for custom calls to `ARKStepSet...` methods.
     */
    std::function<void(void *arkode_mem)> custom_setup;

  private:
    /**
     * Set up the (non)linear solver and preconditioners in the ARKODE memory
     * object based on the user-specified functions.
     * @param solution The solution vector which is used as a template to create
     *   new vectors.
     */
    void
    setup_system_solver(const VectorType           &solution,
                        internal::InvocationContext inv_ctx);

    /**
     * Set up the solver and preconditioner for a non-identity @ref GlossMassMatrix "mass matrix" in
     * the ARKODE memory object based on the user-specified functions.
     * @param solution The solution vector which is used as a template to create
     *   new vectors.
     */
    void
    setup_mass_solver(const VectorType           &solution,
                      internal::InvocationContext inv_ctx);

    /**
     * ARKode memory object.
     */
    void *arkode_mem;

    /**
     * ARKStepper configuration data.
     */
    AdditionalData data;

    std::unique_ptr<internal::LinearSolverWrapper<VectorType>> linear_solver;
    std::unique_ptr<internal::LinearSolverWrapper<VectorType>> mass_solver;

    internal::ARKCallbackContext<ARKStepper<VectorType>> callback_ctx;
  };


  template <typename VectorType>
  ARKStepper<VectorType>::AdditionalData::AdditionalData(
    const unsigned int maximum_non_linear_iterations,
    const bool         implicit_function_is_linear,
    const bool         implicit_function_is_time_independent,
    const bool         mass_is_time_independent,
    const int          anderson_acceleration_subspace)
    : maximum_non_linear_iterations(maximum_non_linear_iterations)
    , implicit_function_is_linear(implicit_function_is_linear)
    , implicit_function_is_time_independent(
        implicit_function_is_time_independent)
    , mass_is_time_independent(mass_is_time_independent)
    , anderson_acceleration_subspace(anderson_acceleration_subspace)
  {}



  template <typename VectorType>
  void
  ARKStepper<VectorType>::AdditionalData::add_parameters(ParameterHandler &prm)
  {
    prm.add_parameter("Maximum number of nonlinear iterations",
                      maximum_non_linear_iterations);
    prm.add_parameter("Implicit function is linear",
                      implicit_function_is_linear);
    prm.add_parameter("Implicit function is time independent",
                      implicit_function_is_time_independent);
    prm.add_parameter("Mass is time independent", mass_is_time_independent);
    prm.add_parameter("Anderson-acceleration subspace",
                      anderson_acceleration_subspace);
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
