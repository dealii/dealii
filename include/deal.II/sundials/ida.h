//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2020 by the deal.II authors
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

#ifndef dealii_sundials_ida_h
#define dealii_sundials_ida_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#ifdef DEAL_II_WITH_SUNDIALS

#  include <deal.II/base/conditional_ostream.h>
#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/logstream.h>
#  include <deal.II/base/parameter_handler.h>
#  ifdef DEAL_II_WITH_PETSC
#    include <deal.II/lac/petsc_block_vector.h>
#    include <deal.II/lac/petsc_vector.h>
#  endif
#  include <deal.II/lac/vector.h>
#  include <deal.II/lac/vector_memory.h>

#  ifdef DEAL_II_SUNDIALS_WITH_IDAS
#    include <idas/idas.h>
#  else
#    include <ida/ida.h>
#  endif

#  include <sundials/sundials_config.h>
#  if DEAL_II_SUNDIALS_VERSION_LT(3, 0, 0)
#    include <ida/ida_spbcgs.h>
#    include <ida/ida_spgmr.h>
#    include <ida/ida_sptfqmr.h>
#  endif
#  include <boost/signals2.hpp>

#  include <nvector/nvector_serial.h>
#  include <sundials/sundials_math.h>
#  include <sundials/sundials_types.h>

#  include <memory>


DEAL_II_NAMESPACE_OPEN

// Shorthand notation for IDA error codes.
#  define AssertIDA(code) Assert(code >= 0, ExcIDAError(code))

namespace SUNDIALS
{
  /**
   * Interface to SUNDIALS Implicit Differential-Algebraic (IDA) solver.
   *
   * The class IDA is a wrapper to SUNDIALS Implicit Differential-Algebraic
   * solver which is a general purpose solver for systems of
   * Differential-Algebraic Equations (DAEs).
   *
   * The user has to provide the implementation of the following std::functions:
   *  - reinit_vector;
   *  - residual;
   *  - setup_jacobian;
   *  - solve_jacobian_system;
   *
   * Optionally, also the following functions could be rewritten. By default
   * they do nothing, or are not required. If you call the constructor in a way
   * that requires a not-implemented function, an Assertion will be
   * thrown.
   *  - solver_should_restart;
   *  - differential_components;
   *  - get_local_tolerances;
   *
   * To output steps, connect a function to the signal
   *  - output_step;
   *
   * Citing from the SUNDIALS documentation:
   *
   *   Consider a system of Differential-Algebraic Equations written in the
   *   general form
   *
   * \f[
   *   \begin{cases}
   *       F(t,y,\dot y) = 0\, , \\
   *       y(t_0) = y_0\, , \\
   *       \dot y (t_0) = \dot y_0\, .
   *   \end{cases}
   * \f]
   *
   * where $y,\dot y$ are vectors in $\mathbb{R}^n$, $t$ is often the time (but
   * can also be a parametric quantity), and
   * $F:\mathbb{R}\times\mathbb{R}^n\times \mathbb{R}^n\rightarrow\mathbb{R}^n$.
   * Such problem is solved using Newton iteration augmented with a line search
   * global strategy. The integration method used in IDA is the variable-order,
   * variable-coefficient BDF (Backward Differentiation Formula), in
   * fixed-leading-coefficient. The method order ranges from 1 to 5, with
   * the BDF of order $q$ given by the multistep formula
   *
   * \f[
   *   \sum_{i=0}^q \alpha_{n,i}\,y_{n-i}=h_n\,\dot y_n\, ,
   *   \label{eq:bdf}
   * \f]
   *
   * where $y_n$ and $\dot y_n$ are the computed approximations of $y(t_n)$
   * and $\dot y(t_n)$, respectively, and the step size is
   * $h_n=t_n-t_{n-1}$. The coefficients $\alpha_{n,i}$ are uniquely
   * determined by the order $q$, and the history of the step sizes. The
   * application of the BDF method to the DAE system results in a nonlinear
   *algebraic system to be solved at each time step:
   *
   * \f[
   *   G(y_n)\equiv F\left(t_n,y_n,\dfrac{1}{h_n}\sum_{i=0}^q
   *  \alpha_{n,i}\,y_{n-i}\right)=0\, .
   * \f]
   * The Newton method leads to a linear system of the form
   * \f[
   *   J[y_{n(m+1)}-y_{n(m)}]=-G(y_{n(m)})\, ,
   * \f]
   *
   *where $y_{n(m)}$ is the $m$-th approximation to $y_n$, and $J$ is the
   *approximation of the system Jacobian
   *
   * \f[
   *   J=\dfrac{\partial G}{\partial y} = \dfrac{\partial F}{\partial y} +
   *  \alpha \dfrac{\partial F}{\partial \dot y}\, ,
   * \f]
   *
   * and $\alpha = \alpha_{n,0}/h_n$. It is worth mentioning that the
   * scalar $\alpha$ changes whenever the step size or method order
   * changes.
   *
   * To provide a simple example, consider the following harmonic oscillator
   *problem: \f[ \begin{split}
   *   u'' & = -k^2 u \\
   *   u (0) & = 0 \\
   *   u'(0) & = k
   * \end{split}
   * \f]
   *
   * We write it in terms of a first order ode:
   *\f[
   * \begin{matrix}
   *   y_0' & -y_1      & = 0 \\
   *   y_1' & + k^2 y_0 & = 0
   * \end{matrix}
   * \f]
   *
   * That is $F(y', y, t) = y' + A y = 0 $
   * where
   * \f[
   * \begin{matrix}
   * 0 & -1 \\
   * k^2 &0
   * \end{matrix}
   * \f]
   * and $y(0)=(0, k)$, $y'(0) = (k, 0)$.
   *
   * The exact solution is $y_0(t) = \sin(k t)$, $y_1(t) = y_0'(t) = k \cos(k
   *t)$, $y_1'(t) = -k^2 \sin(k t)$.
   *
   * The Jacobian to assemble is the following:  $J = \alpha I + A$.
   *
   * This is achieved by the following snippet of code:
   * @code
   * using VectorType = Vector<double>;
   *
   * VectorType y(2);
   * VectorType y_dot(2);
   *
   * double kappa = 1.0;
   *
   * FullMatrix<double> A(2,2);
   * A(0,1) = -1.0;
   * A(1,0) = kappa*kappa;
   *
   * FullMatrix<double> J(2,2);
   * FullMatrix<double> Jinv(2,2);
   *
   * IDA time_stepper;
   *
   * time_stepper.reinit_vector = [&] (VectorType&v)
   * {
   *   v.reinit(2);
   * };
   *
   * time_stepper.residual = [&](const double t,
   *                             const VectorType &y,
   *                             const VectorType &y_dot,
   *                             VectorType &res) ->int
   * {
   *   res = y_dot;
   *   A.vmult_add(res, y);
   *   return 0;
   * };
   *
   * time_stepper.setup_jacobian = [&](const double ,
   *                                   const VectorType &,
   *                                   const VectorType &,
   *                                   const double alpha) ->int
   * {
   *   J = A;
   *
   *   J(0,0) = alpha;
   *   J(1,1) = alpha;
   *
   *   Jinv.invert(J);
   *   return 0;
   * };
   *
   * time_stepper.solve_jacobian_system = [&](const VectorType &src,
   *                                          VectorType &dst) ->int
   * {
   *   Jinv.vmult(dst,src);
   *   return 0;
   * };
   *
   * y[1] = kappa;
   * y_dot[0] = kappa;
   * time_stepper.solve_dae(y,y_dot);
   * @endcode
   *
   * @author Luca Heltai, Alberto Sartori, 2017.
   */
  template <typename VectorType = Vector<double>>
  class IDA
  {
  public:
    /**
     * Additional parameters that can be passed to the IDA class.
     */
    class AdditionalData
    {
    public:
      /**
       * IDA is a Differential Algebraic solver. As such, it requires initial
       * conditions also for the first order derivatives. If you do not provide
       * consistent initial conditions, (i.e., conditions for which F(y_dot(0),
       * y(0), 0) = 0), you can ask SUNDIALS to compute initial conditions for
       * you by specifying InitialConditionCorrection for the initial
       * conditions both at the `initial_time` (`ic_type`) and after a reset
       * has occurred (`reset_type`).
       */
      enum InitialConditionCorrection
      {
        /**
         * Do not try to make initial conditions consistent.
         */
        none = 0,

        /**
         * Compute the algebraic components of y and differential
         * components of y_dot, given the differential components of y.
         *    This option requires that the user specifies differential and
         *    algebraic components in the function get_differential_components.
         */
        use_y_diff = 1,

        /**
         * Compute all components of y, given y_dot.
         */
        use_y_dot = 2
      };

      /**
       * Initialization parameters for IDA.
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
       * @param maximum_order Maximum BDF order
       * @param maximum_non_linear_iterations Maximum number of nonlinear
       * iterations
       *
       * Error parameters:
       *
       * @param absolute_tolerance Absolute error tolerance
       * @param relative_tolerance Relative error tolerance
       * @param ignore_algebraic_terms_for_errors Ignore algebraic terms for
       * error computations
       *
       * Initial condition correction parameters:
       *
       * @param ic_type Initial condition correction type
       * @param reset_type Initial condition correction type after restart
       * @param maximum_non_linear_iterations_ic Initial condition Newton max
       * iterations
       */
      AdditionalData( // Initial parameters
        const double initial_time      = 0.0,
        const double final_time        = 1.0,
        const double initial_step_size = 1e-2,
        const double output_period     = 1e-1,
        // Running parameters
        const double       minimum_step_size             = 1e-6,
        const unsigned int maximum_order                 = 5,
        const unsigned int maximum_non_linear_iterations = 10,
        // Error parameters
        const double absolute_tolerance                = 1e-6,
        const double relative_tolerance                = 1e-5,
        const bool   ignore_algebraic_terms_for_errors = true,
        // Initial conditions parameters
        const InitialConditionCorrection &ic_type    = use_y_diff,
        const InitialConditionCorrection &reset_type = use_y_diff,
        const unsigned int                maximum_non_linear_iterations_ic = 5)
        : initial_time(initial_time)
        , final_time(final_time)
        , initial_step_size(initial_step_size)
        , minimum_step_size(minimum_step_size)
        , absolute_tolerance(absolute_tolerance)
        , relative_tolerance(relative_tolerance)
        , maximum_order(maximum_order)
        , output_period(output_period)
        , ignore_algebraic_terms_for_errors(ignore_algebraic_terms_for_errors)
        , ic_type(ic_type)
        , reset_type(reset_type)
        , maximum_non_linear_iterations_ic(maximum_non_linear_iterations_ic)
        , maximum_non_linear_iterations(maximum_non_linear_iterations)
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
       *   set Maximum order of BDF                   = 5
       *   set Minimum step size                      = 0.000001
       * end
       * @endcode
       *
       * These are one-to-one with the options you can pass at construction
       * time.
       *
       * The options you pass at construction time are set as default values in
       * the ParameterHandler object `prm`. You can later modify them by parsing
       * a parameter file using `prm`. The values of the parameter will be
       * updated whenever the content of `prm` is updated.
       *
       * Make sure that this class lives longer than `prm`. Undefined behaviour
       * will occur if you destroy this class, and then parse a parameter file
       * using `prm`.
       */
      void
      add_parameters(ParameterHandler &prm)
      {
        prm.add_parameter("Initial time", initial_time);
        prm.add_parameter("Final time", final_time);
        prm.add_parameter("Time interval between each output", output_period);

        prm.enter_subsection("Running parameters");
        prm.add_parameter("Initial step size", initial_step_size);
        prm.add_parameter("Minimum step size", minimum_step_size);
        prm.add_parameter("Maximum order of BDF", maximum_order);
        prm.add_parameter("Maximum number of nonlinear iterations",
                          maximum_non_linear_iterations);
        prm.leave_subsection();

        prm.enter_subsection("Error control");
        prm.add_parameter("Absolute error tolerance", absolute_tolerance);
        prm.add_parameter("Relative error tolerance", relative_tolerance);
        prm.add_parameter(
          "Ignore algebraic terms for error computations",
          ignore_algebraic_terms_for_errors,
          "Indicate whether or not to suppress algebraic variables "
          "in the local error test.");
        prm.leave_subsection();

        prm.enter_subsection("Initial condition correction parameters");
        static std::string ic_type_str = "use_y_diff";
        prm.add_parameter(
          "Correction type at initial time",
          ic_type_str,
          "This is one of the following three options for the "
          "initial condition calculation. \n"
          " none: do not try to make initial conditions consistent. \n"
          " use_y_diff: compute the algebraic components of y and differential\n"
          "    components of y_dot, given the differential components of y. \n"
          "    This option requires that the user specifies differential and \n"
          "    algebraic components in the function get_differential_components.\n"
          " use_y_dot: compute all components of y, given y_dot.",
          Patterns::Selection("none|use_y_diff|use_y_dot"));
        prm.add_action("Correction type at initial time",
                       [&](const std::string &value) {
                         if (value == "use_y_diff")
                           ic_type = use_y_diff;
                         else if (value == "use_y_dot")
                           ic_type = use_y_dot;
                         else if (value == "none")
                           ic_type = none;
                         else
                           AssertThrow(false, ExcInternalError());
                       });

        static std::string reset_type_str = "use_y_diff";
        prm.add_parameter(
          "Correction type after restart",
          reset_type_str,
          "This is one of the following three options for the "
          "initial condition calculation. \n"
          " none: do not try to make initial conditions consistent. \n"
          " use_y_diff: compute the algebraic components of y and differential\n"
          "    components of y_dot, given the differential components of y. \n"
          "    This option requires that the user specifies differential and \n"
          "    algebraic components in the function get_differential_components.\n"
          " use_y_dot: compute all components of y, given y_dot.",
          Patterns::Selection("none|use_y_diff|use_y_dot"));
        prm.add_action("Correction type after restart",
                       [&](const std::string &value) {
                         if (value == "use_y_diff")
                           reset_type = use_y_diff;
                         else if (value == "use_y_dot")
                           reset_type = use_y_dot;
                         else if (value == "none")
                           reset_type = none;
                         else
                           AssertThrow(false, ExcInternalError());
                       });
        prm.add_parameter("Maximum number of nonlinear iterations",
                          maximum_non_linear_iterations_ic);
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
       * Maximum order of BDF.
       */
      unsigned int maximum_order;

      /**
       * Time period between each output.
       */
      double output_period;

      /**
       * Ignore algebraic terms for errors.
       */
      bool ignore_algebraic_terms_for_errors;

      /**
       * Type of correction for initial conditions.
       *
       * If you do not provide consistent initial conditions, (i.e., conditions
       * for which $F(y_dot(0), y(0), 0) = 0$), you can ask SUNDIALS to compute
       * initial conditions for you by using the `ic_type` parameter at
       * construction time.
       *
       * Notice that you could in principle use this capabilities to solve for
       * steady state problems by setting y_dot to zero, and asking to compute
       * $y(0)$ that satisfies $F(0, y(0), 0) = 0$, however the nonlinear solver
       * used inside IDA may not be robust enough for complex problems with
       * several millions unknowns.
       */
      InitialConditionCorrection ic_type;

      /**
       * Type of correction for initial conditions to be used after a solver
       * restart.
       *
       * If you do not have consistent initial conditions after a restart,
       * (i.e., conditions for which F(y_dot(t_restart), y(t_restart),
       * t_restart) = 0), you can ask SUNDIALS to compute the new initial
       * conditions for you by using the `reset_type` parameter at construction
       * time.
       */
      InitialConditionCorrection reset_type;

      /**
       * Maximum number of iterations for Newton method in IC calculation.
       */
      unsigned maximum_non_linear_iterations_ic;

      /**
       * Maximum number of iterations for Newton method during time advancement.
       */
      unsigned int maximum_non_linear_iterations;
    };

    /**
     * Constructor. It is possible to fine tune the SUNDIALS IDA solver by
     * passing an AdditionalData() object that sets all of the solver
     * parameters.
     *
     * IDA is a Differential Algebraic solver. As such, it requires initial
     * conditions also for the first order derivatives. If you do not provide
     * consistent initial conditions, (i.e., conditions for which F(y_dot(0),
     * y(0), 0) = 0), you can ask SUNDIALS to compute initial conditions for you
     * by using the `ic_type` parameter at construction time.
     *
     * You have three options
     * -  none: do not try to make initial conditions consistent.
     * -  use_y_diff: compute the algebraic components of y and differential
     *    components of y_dot, given the differential components of y.
     *    This option requires that the user specifies differential and
     *    algebraic components in the function get_differential_components.
     * -  use_y_dot: compute all components of y, given y_dot.
     *
     * By default, this class assumes that all components are differential, and
     * that you want to solve a standard ode. In this case, the initial
     * component type is set to `use_y_diff`, so that the `y_dot` at time
     * t=`initial_time` is computed by solving the nonlinear problem $F(y_dot,
     * y(t0), t0) = 0$ in the variable `y_dot`.
     *
     * Notice that a Newton solver is used for this computation. The Newton
     * solver parameters can be tweaked by acting on `ic_alpha` and
     * `ic_max_iter`.
     *
     * If you reset the solver at some point, you may want to select a different
     * computation for the initial conditions after reset. Say, for example,
     * that you have refined a grid, and after transferring the solution to the
     * new grid, the initial conditions are no longer consistent. Then you can
     * choose how these are made consistent, using the same three options that
     * you used for the initial conditions in `reset_type`.
     *
     * The MPI communicator is simply ignored in the serial case.
     *
     * @param data IDA configuration data
     * @param mpi_comm MPI communicator
     */
    IDA(const AdditionalData &data     = AdditionalData(),
        const MPI_Comm        mpi_comm = MPI_COMM_WORLD);

    /**
     * Destructor.
     */
    ~IDA();

    /**
     * Integrate differential-algebraic equations. This function returns the
     * final number of computed steps.
     */
    unsigned int
    solve_dae(VectorType &solution, VectorType &solution_dot);

    /**
     * Clear internal memory and start with clean objects. This function is
     * called when the simulation start and when the user returns true to a
     * call to solver_should_restart().
     *
     * By default solver_should_restart() returns false. If the user needs to
     * implement, for example, local adaptivity in space, he or she may assign
     * a different function to solver_should_restart() that performs all mesh
     * changes, transfers the solution and the solution dot to the new mesh,
     * and returns true.
     *
     * During reset(), both y and yp are checked for consistency, and according
     * to what was specified as ic_type (if t==initial_time) or reset_type (if
     * t>initial_time), yp, y, or both are modified to obtain a consistent set
     * of initial data.
     *
     * @param[in] t  The new starting time
     * @param[in] h  The new (tentative) starting time step
     * @param[in,out] y   The new (tentative) initial solution
     * @param[in,out] yp  The new (tentative) initial solution_dot
     */
    void
    reset(const double t, const double h, VectorType &y, VectorType &yp);

    /**
     * Reinit vector to have the right size, MPI communicator, etc.
     */
    std::function<void(VectorType &)> reinit_vector;

    /**
     * Compute residual. Return $F(t, y, \dot y)$.
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (IDAReinit will be called if this happens, and
     *       then last function will be attempted again
     * - <0: Unrecoverable error the computation will be aborted and an
     * assertion will be thrown.
     */
    std::function<int(const double      t,
                      const VectorType &y,
                      const VectorType &y_dot,
                      VectorType &      res)>
      residual;

    /**
     * Compute Jacobian. This function is called by IDA any time a Jacobian
     * update is required. The user should compute the Jacobian (or update all
     * the variables that allow the application of the Jacobian). This function
     * is called by IDA once, before any call to solve_jacobian_system().
     *
     * The Jacobian $J$ should be a (possibly inexact) computation of
     * \f[
     *   J=\dfrac{\partial G}{\partial y} = \dfrac{\partial F}{\partial y} +
     *  \alpha \dfrac{\partial F}{\partial \dot y}.
     * \f]
     *
     * If the user uses a matrix based computation of the Jacobian, than this
     * is the right place where an assembly routine should be called to
     * assemble both a matrix and a preconditioner for the Jacobian system.
     * Subsequent calls (possibly more than one) to solve_jacobian_system() can
     * assume that this function has been called at least once.
     *
     * Notice that no assumption is made by this interface on what the user
     * should do in this function. IDA only assumes that after a call to
     * setup_jacobian() it is possible to call solve_jacobian_system(), to
     * obtain a solution $x$ to the system $J x = b$.
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (IDAReinit will be called if this happens, and
     *       then last function will be attempted again
     * - <0: Unrecoverable error the computation will be aborted and an
     * assertion will be thrown.
     */
    std::function<int(const double      t,
                      const VectorType &y,
                      const VectorType &y_dot,
                      const double      alpha)>
      setup_jacobian;

    /**
     * Solve the Jacobian linear system. This function will be called by IDA
     * (possibly several times) after setup_jacobian() has been called at least
     * once. IDA tries to do its best to call setup_jacobian() the minimum
     * amount of times. If convergence can be achieved without updating the
     * Jacobian, then IDA does not call setup_jacobian() again. If, on the
     * contrary, internal IDA convergence tests fail, then IDA calls again
     * setup_jacobian() with updated vectors and coefficients so that successive
     * calls to solve_jacobian_systems() lead to better convergence in the
     * Newton process.
     *
     * The jacobian $J$ should be (an approximation of) the system Jacobian
     * \f[
     *   J=\dfrac{\partial G}{\partial y} = \dfrac{\partial F}{\partial y} +
     *  \alpha \dfrac{\partial F}{\partial \dot y}.
     * \f]
     *
     * A call to this function should store in `dst` the result of $J^{-1}$
     * applied to `src`, i.e., `J*dst = src`. It is the users responsibility
     * to set up proper solvers and preconditioners inside this function.
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (IDAReinit will be called if this happens, and
     *       then last function will be attempted again
     * - <0: Unrecoverable error the computation will be aborted and an
     * assertion will be thrown.
     */
    std::function<int(const VectorType &rhs, VectorType &dst)>
      solve_jacobian_system;

    /**
     * Process solution. This function is called by IDA at fixed time steps,
     * every `output_period` seconds, and it is passed a polynomial
     * interpolation of the solution and of its time derivative, computed using
     * the current BDF order and the (internally stored) previously computed
     * solution steps.
     *
     * Notice that it is well possible that internally IDA computes a time step
     * which is much larger than the `output_period` step, and therefore calls
     * this function consecutively several times by simply performing all
     * intermediate interpolations. There is no relationship between how many
     * times this function is called and how many time steps have actually been
     * computed.
     */
    std::function<void(const double       t,
                       const VectorType & sol,
                       const VectorType & sol_dot,
                       const unsigned int step_number)>
      output_step;

    /**
     * Evaluate whether the solver should be restarted (for example because the
     * number of degrees of freedom has changed).
     *
     * This function is supposed to perform all operations that are necessary in
     * `sol` and `sol_dot` to make sure that the resulting vectors are
     * consistent, and of the correct final size.
     *
     * For example, one may decide that a local refinement is necessary at time
     * t. This function should then return true, and change the dimension of
     * both sol and sol_dot to reflect the new dimension. Since IDA does not
     * know about the new dimension, an internal reset is necessary.
     *
     * The default implementation simply returns `false`, i.e., no restart is
     * performed during the evolution.
     */
    std::function<bool(const double t, VectorType &sol, VectorType &sol_dot)>
      solver_should_restart;

    /**
     * Return an index set containing the differential components.
     * Implementation of this function is optional. The default is to return a
     * complete index set. If your equation is also algebraic (i.e., it
     * contains algebraic constraints, or Lagrange multipliers), you should
     * overwrite this function in order to return only the differential
     * components of your system.
     */
    std::function<IndexSet()> differential_components;

    /**
     * Return a vector whose components are the weights used by IDA to compute
     * the vector norm. The implementation of this function is optional. If the
     * user does not provide an implementation, the weights are assumed to be
     * all ones.
     */
    std::function<VectorType &()> get_local_tolerances;

    /**
     * Handle IDA exceptions.
     */
    DeclException1(ExcIDAError,
                   int,
                   << "One of the SUNDIALS IDA internal functions "
                   << " returned a negative error code: " << arg1
                   << ". Please consult SUNDIALS manual.");


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
     * This function is executed at construction time to set the
     * std::function above to trigger an assert if they are not
     * implemented.
     */
    void
    set_functions_to_trigger_an_assert();

    /**
     * IDA configuration data.
     */
    AdditionalData data;

    /**
     * IDA memory object.
     */
    void *ida_mem;

    /**
     * IDA solution vector.
     */
    N_Vector yy;

    /**
     * IDA solution derivative vector.
     */
    N_Vector yp;

    /**
     * IDA absolute tolerances vector.
     */
    N_Vector abs_tolls;

    /**
     * IDA differential components vector.
     */
    N_Vector diff_id;

    /**
     * MPI communicator. SUNDIALS solver runs happily in
     * parallel. Note that if the library is compiled without MPI
     * support, MPI_Comm is aliased as int.
     */
    MPI_Comm communicator;

    /**
     * Memory pool of vectors.
     */
    GrowingVectorMemory<VectorType> mem;

#  ifdef DEAL_II_WITH_PETSC
#    ifdef PETSC_USE_COMPLEX
    static_assert(!std::is_same<VectorType, PETScWrappers::MPI::Vector>::value,
                  "Sundials does not support complex scalar types, "
                  "but PETSc is configured to use a complex scalar type!");

    static_assert(
      !std::is_same<VectorType, PETScWrappers::MPI::BlockVector>::value,
      "Sundials does not support complex scalar types, "
      "but PETSc is configured to use a complex scalar type!");
#    endif // PETSC_USE_COMPLEX
#  endif   // DEAL_II_WITH_PETSC
  };

} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SUNDIALS

#endif
