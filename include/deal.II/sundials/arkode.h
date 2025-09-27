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
  /**
   * Interface to SUNDIALS additive Runge-Kutta methods (ARKode).
   *
   * The class ARKode is a wrapper to SUNDIALS variable-step, embedded,
   * additive Runge-Kutta solver which is a general purpose solver for systems
   * of ordinary differential equations characterized by the presence of both
   * fast and slow dynamics.
   *
   * Fast dynamics are treated implicitly, and slow dynamics are treated
   * explicitly, using nested families of implicit and explicit Runge-Kutta
   * solvers.
   *
   * Citing directly from ARKode documentation:
   *
   * ARKode solves ODE initial value problems (IVPs) in $R^N$. These problems
   * should be posed in explicit form as
   *
   * \f[
   *   M\dot y = f_E(t, y) + f_I (t, y), \qquad y(t_0) = y_0.
   * \f]
   *
   * Here, $t$ is the independent variable (e.g. time), and the dependent
   * variables are given by $y \in R^N$, and we use notation $\dot y$ to denote
   * $dy/dt$. $M$ is a user-supplied nonsingular operator from $R^N \to R^N$.
   * This operator may depend on $t$ but not on $y$.
   *
   * For standard systems of ordinary differential equations and for problems
   * arising from the spatial semi-discretization of partial differential
   * equations using finite difference or finite volume methods, $M$ is
   * typically the identity matrix, $I$. For PDEs using a finite-element
   * spatial semi-discretization $M$ is typically a well-conditioned mass
   * matrix.
   *
   * The two right-hand side functions may be described as:
   *
   * - $f_E(t, y)$: contains the "slow" time scale components of the system.
   *                This will be integrated using explicit methods.
   * - $f_I(t, y)$: contains the "fast" time scale components of the system.
   *                This will be integrated using implicit methods.
   *
   * ARKode may be used to solve stiff, nonstiff and multi-rate problems.
   * Roughly speaking, stiffness is characterized by the presence of at least
   * one rapidly damped mode, whose time constant is small compared to the time
   * scale of the solution itself. In the implicit/explicit (ImEx) splitting
   * above, these stiff components should be included in the right-hand side
   * function $f_I (t, y)$.
   *
   * For multi-rate problems, a user should provide both of the functions $f_E$
   * and $f_I$ that define the IVP system.
   *
   * For nonstiff problems, only $f_E$ should be provided, and $f_I$ is assumed
   * to be zero, i.e. the system reduces to the non-split IVP:
   *
   * \f[
   *   M\dot y = f_E(t, y), \qquad y(t_0) = y_0.
   * \f]
   *
   * In this scenario, the ARK methods reduce to classical explicit Runge-Kutta
   * methods (ERK). For these classes of methods, ARKode allows orders of
   * accuracy $q = \{2, 3, 4, 5, 6, 8\}$, with embeddings of orders $p = \{1,
   * 2, 3, 4, 5, 7\}$. These default to the Heun-Euler-2-1-2,
   * Bogacki-Shampine-4-2-3, Zonneveld-5-3-4, Cash-Karp-6-4-5, Verner-8-5-6 and
   * Fehlberg-13-7-8 methods, respectively.
   *
   * Finally, for stiff (linear or nonlinear) problems the user may provide only
   * $f_I$, implying that $f_E = 0$, so that the system reduces to the non-split
   * IVP
   *
   * \f[
   *   M\dot y = f_I(t, y), \qquad y(t_0) = y_0.
   * \f]
   *
   * Similarly to ERK methods, in this scenario the ARK methods reduce to
   * classical diagonally-implicit Runge-Kutta methods (DIRK). For these
   * classes of methods, ARKode allows orders of accuracy $q = \{2, 3, 4, 5\}$,
   * with embeddings of orders $p = \{1, 2, 3, 4\}$. These default to the
   * SDIRK-2-1-2, ARK-4-2-3 (implicit), SDIRK-5-3-4 and ARK-8-4-5 (implicit)
   * methods, respectively.
   *
   * For both DIRK and ARK methods, an implicit system of the form
   * \f[
   *  G(z_i) \dealcoloneq M z_i - h_n A^I_{i,i} f_I (t^I_{n,i}, z_i) - a_i = 0
   * \f]
   * must be solved for each stage $z_i , i = 1, \ldots, s$, where
   * we have the data
   * \f[
   *  a_i \dealcoloneq
   *  M y_{n-1} + h_n \sum_{j=1}^{i-1} [ A^E_{i,j} f_E(t^E_{n,j}, z_j)
   *  + A^I_{i,j} f_I (t^I_{n,j}, z_j)]
   * \f]
   * for the ARK methods, or
   * \f[
   *  a_i \dealcoloneq
   *  M y_{n-1} + h_n \sum_{j=1}^{i-1} A^I_{i,j} f_I (t^I_{n,j}, z_j)
   * \f]
   * for the DIRK methods. Here $A^I_{i,j}$ and $A^E_{i,j}$ are the Butcher's
   * tables for the chosen solver.
   *
   * If $f_I(t,y)$ depends nonlinearly on $y$ then the systems above correspond
   * to a nonlinear system of equations; if $f_I (t, y)$ depends linearly on
   * $y$ then this is a linear system of equations. By specifying the flag
   * `implicit_function_is_linear`, ARKode takes some shortcuts that allow a
   * faster solution process.
   *
   * For systems of either type, ARKode allows a choice of solution strategy.
   * The default solver choice is a variant of Newton's method,
   * \f[
   *  z_i^{m+1} = z_i^m +\delta^{m+1},
   * \f]
   * where $m$ is the Newton step index, and the Newton update $\delta^{m+1}$
   * requires the solution of the linear Newton system
   * \f[
   *  N(z_i^m) \delta^{m+1} = -G(z_i^m),
   * \f]
   * where
   * \f[
   * N \dealcoloneq M - \gamma J, \quad J
   * \dealcoloneq \frac{\partial f_I}{\partial y},
   * \qquad \gamma\dealcoloneq h_n A^I_{i,i}.
   * \f]
   *
   * As an alternate to Newton's method, ARKode may solve for each stage $z_i ,i
   * = 1, \ldots , s$ using an Anderson-accelerated fixed point iteration
   * \f[
   * z_i^{m+1} = g(z_i^{m}), m=0,1,\ldots.
   * \f]
   *
   * Unlike with Newton's method, this option does not require the solution of
   * a linear system at each iteration, instead opting for solution of a
   * low-dimensional least-squares solution to construct the nonlinear update.
   *
   * Finally, if the user specifies `implicit_function_is_linear`, i.e.,
   * $f_I(t, y)$ depends linearly on $y$, and if the Newton-based nonlinear
   * solver is chosen, then the system will be solved using only a single
   * Newton iteration. Notice that in order for the Newton solver to be used,
   * then jacobian_times_vector() should be supplied. If it is not supplied then
   * only the fixed-point iteration will be supported, and the
   * `implicit_function_is_linear` setting is ignored.
   *
   * The optimal solver (Newton vs fixed-point) is highly problem-dependent.
   * Since fixed-point solvers do not require the solution of any linear
   * systems, each iteration may be significantly less costly than their Newton
   * counterparts. However, this can come at the cost of slower convergence (or
   * even divergence) in comparison with Newton-like methods. These fixed-point
   * solvers do allow for user specification of the Anderson-accelerated
   * subspace size, $m_k$. While the required amount of solver memory grows
   * proportionately to $m_k N$, larger values of $m_k$ may result in faster
   * convergence.
   *
   * This improvement may be significant even for "small" values, e.g. $1 \leq
   * m_k \leq 5$, and convergence may not improve (or even deteriorate) for
   * larger values of $m_k$. While ARKode uses a Newton-based iteration as its
   * default solver due to its increased robustness on very stiff problems, it
   * is highly recommended that users also consider the fixed-point solver for
   * their cases when attempting a new problem.
   *
   * For either the Newton or fixed-point solvers, it is well-known that both
   * the efficiency and robustness of the algorithm intimately depends on the
   * choice of a good initial guess. In ARKode, the initial guess for either
   * nonlinear solution method is a predicted value $z_i(0)$ that is computed
   * explicitly from the previously-computed data (e.g. $y_{n-2}, y_{n-1}$, and
   * $z_j$ where $j < i$). Additional information on the specific predictor
   * algorithms implemented in ARKode is provided in ARKode documentation.
   *
   * The user has to provide the implementation of at least one (or both) of the
   * following `std::function`s:
   *  - implicit_function()
   *  - explicit_function()
   *
   * If the @ref GlossMassMatrix "mass matrix" is different from the identity, the user should supply
   *  - mass_times_vector() and, optionally,
   *  - mass_times_setup()
   *
   * If the use of a Newton method is desired, then the user should also supply
   * jacobian_times_vector(). jacobian_times_setup() is optional.
   *
   * @note Although SUNDIALS can provide a difference quotient approximation
   *   of the Jacobian, this is currently not supported through this wrapper.
   *
   * A SUNDIALS default solver (SPGMR) is used to solve the linear systems. To
   * use a custom linear solver for the mass matrix and/or Jacobian, set:
   *  - solve_mass() and/or
   *  - solve_linearized_system()
   *
   * To use a custom preconditioner with either a default or custom linear
   * solver, set:
   * - jacobian_preconditioner_solve() and/or mass_preconditioner_solve()
   * and, optionally,
   * - jacobian_preconditioner_setup() and/or mass_preconditioner_setup()
   *
   * Also the following functions could be rewritten. By default
   * they do nothing, or are not required.
   *  - solver_should_restart()
   *  - get_local_tolerances()
   *
   * To produce output at fixed steps, set the function
   *  - output_step()
   *
   * Any other custom settings of the ARKODE object can be specified in
   *  - custom_setup()
   *
   * To provide a simple example, consider the harmonic oscillator problem:
   * \f[
   * \begin{split}
   *   u'' & = -k^2 u \\
   *   u (0) & = 0 \\
   *   u'(0) & = k
   * \end{split}
   * \f]
   *
   * We write it in terms of a first order ode:
   *\f[
   * \begin{matrix}
   *   y_0' & =  y_1 \\
   *   y_1' & = - k^2 y_0
   * \end{matrix}
   * \f]
   *
   * That is $y' = A y$
   * where
   * \f[
   * A \dealcoloneq
   * \begin{pmatrix}
   * 0 & 1 \\
   * -k^2 &0
   * \end{pmatrix}
   * \f]
   * and $y(0)=(0, k)^T$.
   *
   * The exact solution is $y_0(t) = \sin(k t)$, $y_1(t) = y_0'(t) = k \cos(k
   *t)$, $y_1'(t) = -k^2 \sin(k t)$.
   *
   * A minimal implementation, using only explicit RK methods, is given by the
   * following code snippet:
   *
   * @code
   * using VectorType = Vector<double>;
   *
   * SUNDIALS::ARKode<VectorType> ode;
   *
   * const double k = 1.0;
   *
   * ode.explicit_function = [k] (const double / * time * /,
   *                              const VectorType &y,
   *                              VectorType &ydot)
   * {
   *   ydot[0] = y[1];
   *   ydot[1] = -k*k*y[0];
   * };
   *
   * Vector<double> y(2);
   * y[1] = k;
   *
   * ode.solve_ode(y);
   * @endcode
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
     * Constructor.
     *
     * @param data ARKode configuration data
     * @param mpi_comm MPI Communicator over which logging operations are
     * computed. Only used in SUNDIALS 6 and newer.
     */
    ARKode(const AdditionalData &data, const MPI_Comm mpi_comm);

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
     * Set up the (non)linear solver and preconditioners in the ARKODE memory
     * object based on the user-specified functions.
     * @param solution The solution vector which is used as a template to create
     *   new vectors.
     */
    void
    setup_system_solver(const VectorType &solution);

    /**
     * Set up the solver and preconditioner for a non-identity @ref GlossMassMatrix "mass matrix" in
     * the ARKODE memory object based on the user-specified functions.
     * @param solution The solution vector which is used as a template to create
     *   new vectors.
     */
    void
    setup_mass_solver(const VectorType &solution);

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
     * ARKode memory object.
     */
    void *arkode_mem;

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

    std::unique_ptr<internal::LinearSolverWrapper<VectorType>> linear_solver;
    std::unique_ptr<internal::LinearSolverWrapper<VectorType>> mass_solver;

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


  /**
   * Handle ARKode exceptions.
   */
  DeclException1(ExcARKodeError,
                 int,
                 << "One of the SUNDIALS ARKode internal functions "
                 << " returned a negative error code: " << arg1
                 << ". Please consult SUNDIALS manual.");


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
