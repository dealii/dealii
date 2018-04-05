//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2018 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE at
//    the top level of the deal.II distribution.
//
//---------------------------------------------------------------

#ifndef dealii_sundials_arkode_h
#define dealii_sundials_arkode_h

#include <deal.II/base/config.h>
#include <deal.II/base/mpi.h>

#ifdef DEAL_II_WITH_SUNDIALS

#include <deal.II/base/logstream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/lac/petsc_parallel_block_vector.h>
#  include <deal.II/lac/petsc_parallel_vector.h>
#endif
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>


#include <arkode/arkode.h>
#include <arkode/arkode_impl.h>
#include <nvector/nvector_serial.h>
#ifdef DEAL_II_WITH_MPI
#  include <nvector/nvector_parallel.h>
#endif
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include <boost/signals2.hpp>
#include <memory>


DEAL_II_NAMESPACE_OPEN

// Shorthand notation for ARKODE error codes.
#define AssertARKode(code) Assert(code >= 0, ExcARKodeError(code))

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
   *  G(z_i) := M z_i − h_n A^I_{i,i} f_I (t^I_{n,i}, z_i) − a_i = 0
   * \f]
   * must be solved for each stage $z_i , i = 1, \ldot, s$, where
   * we have the data
   * \[
   *  a_i :=  M y_{n−1} + h_n \sum_{j=1}^{i−1} [ A^E_{i,j} f_E(t^E_{n,j}, z_j)
   *  + A^I_{i,j} f_I (t^I_{n,j}, z_j)]
   * \]
   * for the ARK methods, or
   * \[
   *  a_i :=  M y_{n−1} + h_n \sum_{j=1}^{i−1} A^I_{i,j} f_I (t^I_{n,j}, z_j)
   * \]
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
   * The default solver choice is a variant of Newton’s method,
   * \[
   *  z_i^{m+1} = z_i^m +\delta^{m+1},
   * \]
   * where $m$ is the Newton index, and the Newton update $\delta^{m+1}$
   * requires the solution of the linear Newton system
   * \[
   *  N(z_i^m) \delta^{m+1} = -G(z_i^m),
   * \]
   * where
   * \[
   * N := M - \gamma J, \quad J := \frac{\partial f_I}{\partial y},
   * \qquad \gamma:= h_n A^I_{i,i}.
   * \]
   *
   * As an alternate to Newton’s method, ARKode may solve for each stage $z_i ,i
   * = 1, \ldots , s$ using an Anderson-accelerated fixed point iteration
   * \[
   * z_i^{m+1} = g(z_i^{m}), m=0,1,\ldots.
   * \]
   *
   * Unlike with Newton’s method, this option does not require the solution of
   * a linear system at each iteration, instead opting for solution of a
   * low-dimensional least-squares solution to construct the nonlinear update.
   *
   * Finally, if the user specifies `implicit_function_is_linear`, i.e.,
   * $f_I(t, y)$ depends linearly on $y$, and if the Newton-based nonlinear
   * solver is chosen, then the system will be solved using only a single
   * Newton iteration. Notice that in order for the Newton solver to be used,
   * at least the solve_jacobian_system() function should be supplied. If this
   * function is not supplied, then only the fixed-point iteration will be
   * supported, and the `implicit_function_is_linear` setting is ignored.
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
   * explicitly from the previously-computed data (e.g. $y_{n−2}, y_{n−1}$, and
   * $z_j$ where $j < i$). Additional information on the specific predictor
   * algorithms implemented in ARKode is provided in ARKode documentation.
   *
   * The user has to provide the implementation of the following std::function:
   *  - reinit_vector;
   * and either one or both of
   *  - implicit_function;
   *  - explicit_function;
   *
   * If the mass matrix is different from the identity, the user should supply
   *  - solve_mass_system;
   * and, optionally,
   *  - setup_mass;
   *
   * If the use of a Newton method is desired, then the user should also supply
   *  - solve_jacobian_system;
   * and optionally
   *  - setup_jacobian;
   *
   * Also the following functions could be rewritten. By default
   * they do nothing, or are not required.
   *  - solver_should_restart;
   *  - get_local_tolerances;
   *
   * To produce output at fixed steps, overload the function
   *  - output_step;
   *
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
   * A:=
   * \begin{matrix}
   * 0 & 1 \\
   * -k^2 &0
   * \end{matrix}
   * \f]
   * and $y(0)=(0, k)$.
   *
   * The exact solution is $y_0(t) = \sin(k t)$, $y_1(t) = y_0'(t) = k \cos(k t)$,
   * $y_1'(t) = -k^2 \sin(k t)$.
   *
   * A minimal implementation, using only explicit RK methods, is given by the
   * following code snippet:
   *
   * @code
   * typedef Vector<double> VectorType;
   *
   * SUNDIALS::ARKode<VectorType> ode;
   *
   * ode.reinit_vector = [] (VectorType&v)
   * {
   *   v.reinit(2);
   * };
   *
   * double kappa = 1.0;
   *
   * ode.explicit_function = [kappa] (double,
   *                                  const VectorType &y,
   *                                  VectorType &ydot) -> int
   * {
   *   ydot[0] = y[1];
   *   ydot[1] = -kappa*kappa*y[0];
   *   return 0;
   * };
   *
   * Vector<double> y(2);
   * y[1] = kappa;
   *
   * ode.solve_ode(y);
   * @endcode
   *
   * @author Luca Heltai, 2017.
   */
  template<typename VectorType=Vector<double> >
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
       * @param output_period Time interval between each output
       *
       * Running parameters:
       *
       * @param minimum_step_size Minimum step size
       * @param maximum_order Maximum ARK order
       * @param maximum_non_linear_iterations Maximum number of nonlinear iterations
       * @param implicit_function_is_linear Specifies that the implicit portion of the problem is linear
       * @param implicit_function_is_time_independent Specifies that the implicit portion of the problem
       *        is linear and time independent
       *
       * Error parameters:
       *
       * @param absolute_tolerance Absolute error tolerance
       * @param relative_tolerance Relative error tolerance
       */
      AdditionalData(
        // Initial parameters
        const double &initial_time = 0.0,
        const double &final_time = 1.0,
        const double &initial_step_size = 1e-2,
        const double &output_period = 1e-1,
        // Running parameters
        const double &minimum_step_size = 1e-6,
        const unsigned int &maximum_order = 5,
        const unsigned int &maximum_non_linear_iterations = 10,
        const bool implicit_function_is_linear = false,
        const bool implicit_function_is_time_independent = false,
        // Error parameters
        const double &absolute_tolerance = 1e-6,
        const double &relative_tolerance = 1e-5) :
        initial_time(initial_time),
        final_time(final_time),
        initial_step_size(initial_step_size),
        minimum_step_size(minimum_step_size),
        absolute_tolerance(absolute_tolerance),
        relative_tolerance(relative_tolerance),
        maximum_order(maximum_order),
        output_period(output_period),
        maximum_non_linear_iterations(maximum_non_linear_iterations),
        implicit_function_is_linear(implicit_function_is_linear),
        implicit_function_is_time_independent(implicit_function_is_time_independent)
      {}

      /**
       * Add all AdditionalData() parameters to the given ParameterHandler
       * object. When the parameters are parsed from a file, the internal
       * parameters are automatically updated.
       *
       * The following parameters are declared:
       *
       * @code
       * set Final time                        = 1.000000
       * set Initial time                      = 0.000000
       * set Time interval between each output = 0.2
       * subsection Error control
       *   set Absolute error tolerance                      = 0.000001
       *   set Ignore algebraic terms for error computations = true
       *   set Relative error tolerance                      = 0.00001
       *   set Use local tolerances                          = false
       * end
       * subsection Initial condition correction parameters
       *   set Correction type at initial time        = none
       *   set Correction type after restart          = none
       *   set Maximum number of nonlinear iterations = 5
       * end
       * subsection Running parameters
       *   set Initial step size                      = 0.1
       *   set Maximum number of nonlinear iterations = 10
       *   set Maximum order of ARK                   = 5
       *   set Minimum step size                      = 0.000001
       * end
       * @endcode
       *
       * These are one-to-one with the options you can pass at construction time.
       *
       * The options you pass at construction time are set as default values in
       * the ParameterHandler object `prm`. You can later modify them by parsing
       * a parameter file using `prm`. The values of the parameter will be updated
       * whenever the content of `prm` is updated.
       *
       * Make sure that this class lives longer than `prm`. Undefined behaviour
       * will occur if you destroy this class, and then parse a parameter file
       * using `prm`.
       */
      void add_parameters(ParameterHandler &prm)
      {
        prm.add_parameter("Initial time", initial_time);
        prm.add_parameter("Final time", final_time);
        prm.add_parameter("Time interval between each output", output_period);

        prm.enter_subsection("Running parameters");
        prm.add_parameter("Initial step size",initial_step_size);
        prm.add_parameter("Minimum step size", minimum_step_size);
        prm.add_parameter("Maximum order of ARK", maximum_order);
        prm.add_parameter("Maximum number of nonlinear iterations", maximum_non_linear_iterations);
        prm.add_parameter("Implicit function is linear", implicit_function_is_linear);
        prm.add_parameter("Implicit function is time independent", implicit_function_is_time_independent);
        prm.leave_subsection();

        prm.enter_subsection("Error control");
        prm.add_parameter("Absolute error tolerance", absolute_tolerance);
        prm.add_parameter("Relative error tolerance", relative_tolerance);
        prm.leave_subsection();
      }

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
       * Time period between each output.
       */
      double output_period;

      /**
       * Maximum number of iterations for Newton or fixed point method during
       * time advancement.
       */
      unsigned int maximum_non_linear_iterations;

      /**
       * Specifies that the implicit portion of the problem is linear.
       */
      bool implicit_function_is_linear;

      /**
       * Specifies that the implicit portion of the problem is linear and time
       * independent.
       */
      bool implicit_function_is_time_independent;
    };

    /**
     * Constructor. It is possible to fine tune the SUNDIALS ARKode solver by
     * passing an AdditionalData() object that sets all of the solver
     * parameters.
     *
     * The MPI communicator is simply ignored in the serial case.
     *
     *
     * @param data ARKode configuration data
     * @param mpi_comm MPI communicator
     */
    ARKode(const AdditionalData &data=AdditionalData(),
           const MPI_Comm mpi_comm = MPI_COMM_WORLD);

    /**
     * Destructor.
     */
    ~ARKode();

    /**
     * Integrate the initial value problem. This function returns the final
     * number of computed steps.
     */
    unsigned int solve_ode(VectorType &solution);

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
     * @param[in] t  The new starting time
     * @param[in] h  The new starting time step
     * @param[in,out] y   The new initial solution
     */
    void reset(const double &t,
               const double &h,
               const VectorType &y);

    /**
     * A function object that users need to supply and that is intended to
     * reinit the given vector.
     */
    std::function<void(VectorType &)> reinit_vector;

    /**
     * A function object that users may supply and that is intended to compute
     * the explicit part of the IVP right hand side. Sets $explicit_f = f_E(t,
     * y)$.
     *
     * At least one of explicit_function() or implicit_function() must be
     * provided. According to which one is provided, explicit, implicit, or
     * mixed RK methods are used.
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (ARKodeReinit will be called if this happens, and
     *       then last function will be attempted again
     * - <0: Unrecoverable error the computation will be aborted and an assertion
     *       will be thrown.
     */
    std::function<int(const double t,
                      const VectorType &y,
                      VectorType &explicit_f)> explicit_function;

    /**
     * A function object that users may supply and that is intended to compute
     * the implicit part of the IVP right hand side. Sets $implicit_f = f_I(t,
     * y)$.
     *
     * At least one of explicit_function() or implicit_function() must be
     * provided. According to which one is provided, explicit, implicit, or
     * mixed RK methods are used.
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (ARKodeReinit will be called if this happens, and
     *       then last function will be attempted again
     * - <0: Unrecoverable error the computation will be aborted and an assertion
     *       will be thrown.
     */
    std::function<int(const double t,
                      const VectorType &y,
                      VectorType &res)> implicit_function;

    /**
     * A function object that users may supply and that is intended to
     * prepare the linear solver for subsequent calls to
     * solve_jacobian_system().
     *
     * Make sure that after a call to this function, we know how to compute
     * solutions of systems $A x = b$, where $A$ is some approximation to the
     * Newton matrix, $M − \gamma \partial f_I/\partial y$. This function is
     * optional. If the user does not provide it, then solve_jacobian_system()
     * is assumed to also perform the setup internally.
     *
     * The setup_jacobian() function may call a user-supplied function to
     * compute needed data related to the Jacobian matrix. Alternatively, it may
     * choose to retrieve and use stored values of this data. In either case,
     * setup_jacobian() may also preprocess that data as needed for
     * solve_jacobian_system(), which may involve calling a generic function
     * (such as for LU factorization).
     *
     * This data may be intended either for direct use (in a direct linear
     * solver) or for use in a preconditioner (in a preconditioned iterative
     * linear solver). The setup_jacobian() function is not called at every
     * stage solve (or even every time step), but only as frequently as the
     * solver determines that it is appropriate to perform the setup task. In
     * this way, Jacobian-related data generated by setup_jacobian() is
     * expected to be used over a number of time steps.
     *
     * If the user uses a matrix based computation of the Jacobian, then this
     * is the right place where an assembly routine shoulde be called to
     * assemble both a matrix and a preconditioner for the Jacobian system.
     * Subsequent calls (possibly more than one) to solve_jacobian_system() can
     * assume that this function has been called at least once.
     *
     * Notice that no assumption is made by this interface on what the user
     * should do in this function. ARKode only assumes that after a call to
     * setup_jacobian() it is possible to call solve_jacobian_system(), to
     * obtain a solution $x$ to the system $J x = b$. If this function is not
     * provided, then it is never called.
     *
     * Arguments to the function are
     *
     * @param[in] t  the current time
     * @param[in] gamma  the current factor to use in the jacobian computation
     * @param[in] ypred  is the predicted $y$ vector for the current ARKode internal step
     * @param[in] fpred  is the value of the implicit right-hand side at ypred,
     *        $f_I (t_n, ypred)$.
     *
     * @param[in] convfail – an input flag used to indicate any problem that occurred
     *   during the solution of the nonlinear equation on the current time step
     *   for which the linear solver is being used. This flag can be used to help
     *   decide whether the Jacobian data kept by a linear solver needs to be
     *   updated or not. Its possible values are:
     *
     *   - ARK_NO_FAILURES: this value is passed if either this is the first call
     *     for this step, or the local error test failed on the previous attempt at
     *     this step (but the Newton iteration converged).
     *
     *   - ARK_FAIL_BAD_J: this value is passed if (a) the previous Newton
     *     corrector iteration did not converge and the linear solver's setup
     *     function indicated that its Jacobian-related data is not current, or (b)
     *     during the previous Newton corrector iteration, the linear solver's
     *     solve function failed in a recoverable manner and the linear solver's
     *     setup function indicated that its Jacobian-related data is not current.
     *
     *   - ARK_FAIL_OTHER: this value is passed if during the current internal
     *     step try, the previous Newton iteration failed to converge even though
     *     the linear solver was using current Jacobian-related data.
     *
     * @param[out] j_is_current: a boolean to be filled in by setup_jacobian(). The value
     * should be set to `true` if the Jacobian data is current after the call,
     * and should be set to `false` if its Jacobian data is not current. If
     * setup_jacobian() calls for re-evaluation of Jacobian data (based on
     * convfail and ARKode state data), then it should set `j_is_current` to
     * `true` unconditionally, otherwise an infinite loop can result.
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (ARKodeReinit will be called if this happens, and
     *       then last function will be attempted again
     * - <0: Unrecoverable error the computation will be aborted and an assertion
     *       will be thrown.
     */
    std::function<int(const int convfail,
                      const double t,
                      const double gamma,
                      const VectorType &ypred,
                      const VectorType &fpred,
                      bool &j_is_current)> setup_jacobian;

    /**
     * A function object that users may supply and that is intended to solve
     * the Jacobian linear system. This function will be called by ARKode
     * (possibly several times) after setup_jacobian() has been called at least
     * once. ARKode tries to do its best to call setup_jacobian() the minimum
     * amount of times. If convergence can be achieved without updating the
     * Jacobian, then ARKode does not call setup_jacobian() again. If, on the
     * contrary, internal ARKode convergence tests fail, then ARKode calls
     * again setup_jacobian() with updated vectors and coefficients so that
     * successive calls to solve_jacobian_systems() lead to better convergence
     * in the Newton process.
     *
     * If you do not specify a solve_jacobian_system() function, then a fixed
     * point iteration is used instead of a Newton method. Notice that this may
     * not converge, or may converge very slowly.
     *
     * The jacobian $J$ should be (an approximation of) the system Jacobian
     * \f[
     *   J = M - \gamma \frac{\partial f_I}{\partial y}
     * \f]
     * evaluated at `t`, `ycur`. `fcur` is $f_I(t,ycur)$.
     *
     * A call to this function should store in `dst` the result of $J^{-1}$
     * applied to `src`, i.e., `J*dst = src`. It is the users responsibility to
     * set up proper solvers and preconditioners inside this function.
     *
     *
     * Arguments to the function are
     *
     * @param[in] t  the current time
     * @param[in] gamma  the current factor to use in the jacobian computation
     * @param[in] ycur  is the current $y$ vector for the current ARKode internal step
     * @param[in] fcur  is the current value of the implicit right-hand side at ycur,
     *        $f_I (t_n, ypred)$.
     *
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (ARKodeReinit will be called if this happens, and
     *       then last function will be attempted again
     * - <0: Unrecoverable error the computation will be aborted and an assertion
     *       will be thrown.
     */
    std::function<int(const double t,
                      const double gamma,
                      const VectorType &ycur,
                      const VectorType &fcur,
                      const VectorType &rhs,
                      VectorType &dst)> solve_jacobian_system;


    /**
     * A function object that users may supply and that is intended to setup
     * the mass matrix. This function is called by ARKode any time a mass
     * matrix update is required. The user should compute the mass matrix (or
     * update all the variables that allow the application of the mass matrix).
     * This function is called by ARKode once, before any call to
     * solve_mass_system().
     *
     * ARKode supports the case where the mass matrix may depend on time, but
     * not the case where the mass matrix depends on the solution itself.
     *
     * If the user does not provide a solve_mass_matrix() function, then the
     * identity is used. If the setup_mass() function is not provided, then
     * solve_mass_system() should do all the work by itself.
     *
     * If the user uses a matrix based computation of the mass matrix, then
     * this is the right place where an assembly routine shoulde be called to
     * assemble both a matrix and a preconditioner. Subsequent calls (possibly
     * more than one) to solve_mass_system() can assume that this function
     * has been called at least once.
     *
     * Notice that no assumption is made by this interface on what the user
     * should do in this function. ARKode only assumes that after a call to
     * setup_mass() it is possible to call solve_mass_system(), to
     * obtain a solution $x$ to the system $M x = b$.
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (ARKodeReinit will be called if this happens, and
     *       then last function will be attempted again
     * - <0: Unrecoverable error the computation will be aborted and an assertion
     *       will be thrown.
     */
    std::function<int(const double t)> setup_mass;

    /**
     * A function object that users may supply and that is intended to solve
     * the mass matrix linear system. This function will be called by ARKode
     * (possibly several times) after setup_mass() has been called at least
     * once. ARKode tries to do its best to call setup_mass() the minimum
     * amount of times.
     *
     * A call to this function should store in `dst` the result of $M^{-1}$
     * applied to `src`, i.e., `M*dst = src`. It is the users responsibility to
     * set up proper solvers and preconditioners inside this function.
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (ARKodeReinit will be called if this happens, and
     *       then last function will be attempted again
     * - <0: Unrecoverable error the computation will be aborted and an assertion
     *       will be thrown.
     */
    std::function<int(const VectorType &rhs, VectorType &dst)> solve_mass_system;

    /**
     * A function object that users may supply and that is intended to
     * postprocess the solution. This function is called by ARKode at fixed
     * time increments (every `output_period` seconds), and it is passed a
     * polynomial interpolation of the solution, computed using the current ARK
     * order and the (internally stored) previously computed solution steps.
     *
     * Notice that it is well possible that internally ARKode computes a time step
     * which is much larger than the `output_period` step, and therefore calls
     * this function consecutively several times by simply performing all
     * intermediate interpolations. There is no relationship between how many
     * times this function is called and how many time steps have actually been
     * computed.
     */
    std::function<void(const double t,
                       const VectorType &sol,
                       const unsigned int step_number)> output_step;

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
    std::function<bool (const double t,
                        VectorType &sol)> solver_should_restart;

    /**
     * A function object that users may supply and that is intended to return a
     * vector whose components are the weights used by ARKode to compute the
     * vector norm. The implementation of this function is optional, and it is
     * used only if implemented.
     */
    std::function<VectorType&()> get_local_tolerances;

    /**
     * Handle ARKode exceptions.
     */
    DeclException1(ExcARKodeError, int, << "One of the SUNDIALS ARKode internal functions "
                   << " returned a negative error code: "
                   << arg1 << ". Please consult SUNDIALS manual.");


  private:

    /**
     * Throw an exception when a function with the given name is not implemented.
     */
    DeclException1(ExcFunctionNotProvided, std::string,
                   << "Please provide an implementation for the function \"" << arg1 << "\"");

    /**
     * This function is executed at construction time to set the
     * std::function above to trigger an assert if they are not
     * implemented.
     */
    void set_functions_to_trigger_an_assert();

    /**
     * ARKode configuration data.
     */
    AdditionalData data;

    /**
     * ARKode memory object.
     */
    void *arkode_mem;

    /**
     * ARKode solution vector.
     */
    N_Vector yy;

    /**
     * ARKode absolute tolerances vector.
     */
    N_Vector abs_tolls;

    /**
     * MPI communicator. SUNDIALS solver runs happily in
     * parallel. Note that if the library is compiled without MPI
     * support, MPI_Comm is typedefed as int.
     */
    MPI_Comm communicator;

    /**
     * Memory pool of vectors.
     */
    GrowingVectorMemory<VectorType> mem;

#ifdef DEAL_II_WITH_PETSC
#ifdef PETSC_USE_COMPLEX
    static_assert(!std::is_same<VectorType, PETScWrappers::MPI::Vector>::value,
                  "Sundials does not support complex scalar types, "
                  "but PETSc is configured to use a complex scalar type!");

    static_assert(!std::is_same<VectorType, PETScWrappers::MPI::BlockVector>::value,
                  "Sundials does not support complex scalar types, "
                  "but PETSc is configured to use a complex scalar type!");
#endif // PETSC_USE_COMPLEX
#endif // DEAL_II_WITH_PETSC
  };

}

DEAL_II_NAMESPACE_CLOSE
#endif


#endif
