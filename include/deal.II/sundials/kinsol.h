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

#ifndef dealii_sundials_kinsol_h
#define dealii_sundials_kinsol_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS


#  include <deal.II/base/conditional_ostream.h>
#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/logstream.h>
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/parameter_handler.h>

#  include <deal.II/lac/vector.h>
#  include <deal.II/lac/vector_memory.h>

#  include <boost/signals2.hpp>

#  include <kinsol/kinsol.h>
#  if DEAL_II_SUNDIALS_VERSION_LT(4, 1, 0)
#    include <kinsol/kinsol_impl.h>
#  endif
#  include <nvector/nvector_serial.h>
#  include <sundials/sundials_math.h>
#  include <sundials/sundials_types.h>

#  include <memory>


DEAL_II_NAMESPACE_OPEN

// Shorthand notation for KINSOL error codes.
#  define AssertKINSOL(code) Assert(code >= 0, ExcKINSOLError(code))

namespace SUNDIALS
{
  /**
   * Interface to SUNDIALS' nonlinear solver (KINSOL).
   *
   * KINSOL is a solver for nonlinear algebraic systems in residual form $F(u)
   * = 0$ or fixed point form $G(u) = u$, where $u$ is a vector which we will
   * assume to be in ${\mathbb R}^n$ or ${\mathbb C}^n$, but that may also have
   * a block structure and may be distributed in parallel computations; the
   * functions $F$ and $G$ satisfy $F,G:{\mathbb R}^N \to{\mathbb R}^N$ or
   * $F,G:{\mathbb C}^N \to{\mathbb C}^N$. It includes a Newton-Krylov solver
   * as well as Picard and fixed point solvers, both of which can be
   * accelerated with Anderson acceleration. KINSOL is based on the previous
   * Fortran package NKSOL of Brown and Saad. An example of using KINSOL
   * can be found in the step-77 tutorial program.
   *
   * KINSOL's Newton solver employs the inexact Newton method. As this solver
   * is intended mainly for large systems, the user is required to provide
   * their own solver function. If a solver function is not provided, the
   * internal dense solver of KINSOL is used. Be warned that this solver
   * computes the Jacobian approximately, and may be efficient only for small
   * systems.
   *
   * At the highest level, KINSOL implements the following iteration
   * scheme:
   *  - set $u_0$ = an initial guess
   *  - For $n = 0, 1, 2, \ldots$ until convergence do:
   *    - Solve $J(u_n)\delta_n = -F(u_n)$
   *    - Set $u_{n+1} = u_n + \lambda \delta_n, 0 < \lambda \leq 1$
   *    - Test for convergence
   *
   * Here, $u_n$ is the $n$-th iterate to $u$, and $J(u) = \nabla_u F(u)$ is
   * the system Jacobian. At each stage in the iteration process, a scalar
   * multiple of the step $\delta_n$, is added to $u_n$ to produce a new
   * iterate, $u_{n+1}$. A test for convergence is made before the iteration
   * continues.
   *
   * Unless specified otherwise by the user, KINSOL strives to update Jacobian
   * information as infrequently as possible to balance the high costs of
   * matrix operations against other costs. Specifically, these updates occur
   * when:
   * - the problem is initialized,
   * - $\|\lambda \delta_{n-1} \|_{D_u,\infty} \geq 1.5$ (inexact Newton only,
   *   see below for a definition of $\| \cdot \|_{D_u,\infty}$)
   * - a specified number of  nonlinear iterations have passed since the last
   *   update,
   * - the linear solver failed recoverably with outdated Jacobian information,
   * - the global strategy failed with outdated Jacobian information, or
   * - $\|\lambda \delta_{n} \|_{D_u,\infty} \leq $ *tolerance* with outdated
   *   Jacobian information.
   *
   * KINSOL allows changes to the above strategy through optional solver
   * inputs. The user can disable the initial Jacobian information evaluation
   * or change the default value of the number of nonlinear iterations after
   * which a Jacobian information update is enforced.
   *
   * To address the case of ill-conditioned nonlinear systems, KINSOL allows
   * prescribing scaling factors both for the solution vector and for the
   * residual vector. For scaling to be used, the user may supply the function
   * get_solution_scaling(), that returns values $D_u$, which are diagonal
   * elements of the scaling matrix such that $D_u u_n$ has all components
   * roughly the same magnitude when $u_n$ is close to a solution, and
   * get_function_scaling(), that supply values $D_F$, which are diagonal
   * scaling matrix elements such that $D_F F$ has all components roughly the
   * same magnitude when $u_n$ is *not* too close to a solution.
   *
   * When scaling values are provided for the solution vector, these values are
   * automatically incorporated into the calculation of the perturbations used
   * for the default difference quotient approximations for Jacobian
   * information if the user does not supply a Jacobian solver through the
   * solve_jacobian_system() function.
   *
   * Two methods of applying a computed step $\delta_n$ to the previously
   * computed solution vector are implemented. The first and simplest is the
   * standard Newton strategy which applies the update with a constant
   * $\lambda$ always set to 1. The other method is a global strategy, which
   * attempts to use the direction implied by $\delta_n$ in the most efficient
   * way for furthering convergence of the nonlinear problem. This technique is
   * implemented in the second strategy, called Linesearch. This option employs
   * both the $\alpha$ and $\beta$ conditions of the Goldstein-Armijo
   * linesearch algorithm given in *J. E. Dennis and R. B. Schnabel. "Numerical
   * Methods for Unconstrained Optimization and Nonlinear Equations." SIAM,
   * Philadelphia, 1996.*, where $\lambda$ is chosen to guarantee a sufficient
   * decrease in $F$ relative to the step length as well as a minimum step
   * length relative to the initial rate of decrease of $F$. One property of the
   * algorithm is that the full Newton step tends to be taken close to the
   * solution.
   *
   * As a user option, KINSOL permits the application of inequality
   * constraints, $u_i > 0$ and $u_i < 0$, as well as $u_i \geq 0$ and $u_i
   * \leq 0$, where $u_i$ is the $i$-th component of $u$. Any such constraint,
   * or no constraint, may be imposed on each component by providing the
   * optional functions
   * - get_lower_than_zero_constrained_entries()
   * - get_greater_than_zero_constrained_entries()
   * - get_lower_equal_than_zero_constrained_entries()
   * - get_greater_or_equal_than_zero_constrained_entries()
   *
   * KINSOL will reduce step lengths in order to ensure that no constraint is
   * violated. Specifically, if a new Newton iterate will violate a constraint,
   * the maximum step length along the Newton direction that will satisfy all
   * constraints is found, and $\delta_n$ is scaled to take a step of that
   * length.
   *
   * The basic fixed-point iteration scheme implemented in KINSOL is given by:
   * - Set $u_0 =$ an initial guess
   * - For $n = 0, 1, 2, \dots$ until convergence do:
   *   - Set $u_{n+1} = G(u_n)$
   *   - Test for convergence
   *
   * At each stage in the iteration process, function $G$ is applied to the
   * current iterate to produce a new iterate, $u_{n+1}$. A test for
   * convergence is made before the iteration continues.
   *
   * For Picard iteration, as implemented in KINSOL, we consider a special form
   * of the nonlinear function $F$, such that $F(u) = Lu - N(u)$, where $L$ is
   * a constant nonsingular matrix and $N$ is (in general) nonlinear.
   *
   * Then the fixed-point function $G$ is defined as $G(u) = u - L^{-1}F(u)$.
   * Within each iteration, the Picard step is computed then added to $u_n$ to
   * produce the new iterate. Next, the nonlinear residual function is
   * evaluated at the new iterate, and convergence is checked. The Picard and
   * fixed point methods can be significantly accelerated using Anderson's
   * method.
   *
   * The user has to provide the implementation of the following std::functions:
   *  - reinit_vector;
   * and only one of
   *  - residual;
   * or
   *  - iteration_function;
   *
   * Specifying residual() allows the user to use Newton strategies (i.e.,
   * $F(u)=0$ will be solved), while specifying iteration_function(), fixed
   * point iteration or Picard iteration will be used (i.e., $G(u)=u$ will be
   * solved).
   *
   * If the use of a Newton method is desired, then the user should also supply
   *  - solve_jacobian_system;
   * and optionally
   *  - setup_jacobian;
   *
   * If the solve_jacobian_system() function is not supplied, then KINSOL will
   * use its internal dense solver for Newton methods, with approximate
   * Jacobian. This may be very expensive for large systems. Fixed point
   * iteration does not require the solution of any linear system.
   *
   * Also the following functions could be rewritten, to provide additional
   * scaling factors for both the solution and the residual evaluation during
   * convergence checks:
   *  - get_solution_scaling;
   *  - get_function_scaling;
   */
  template <typename VectorType = Vector<double>>
  class KINSOL
  {
  public:
    /**
     * Additional parameters that can be passed to the KINSOL class.
     */
    class AdditionalData
    {
    public:
      /**
       * KINSOL solution strategy. KINSOL includes a Newton-Krylov solver (both
       * local and global) as well as Picard and fixed point solvers.
       */
      enum SolutionStrategy
      {
        /**
         * Standard Newton iteration.
         */
        newton = KIN_NONE,
        /**
         * Newton iteration with linesearch.
         */
        linesearch = KIN_LINESEARCH,
        /**
         * Fixed point iteration.
         */
        fixed_point = KIN_FP,
        /**
         * Picard iteration.
         */
        picard = KIN_PICARD,
      };

      /**
       * Initialization parameters for KINSOL.
       *
       * Global parameters:
       *
       * @param strategy Solution strategy
       * @param maximum_non_linear_iterations Maximum number of nonlinear
       * iterations
       * @param function_tolerance %Function norm stopping tolerance
       * @param step_tolerance Scaled step stopping tolerance
       *
       * Newton parameters:
       *
       * @param no_init_setup No initial matrix setup
       * @param maximum_setup_calls Maximum iterations without matrix setup
       * @param maximum_newton_step Maximum allowable scaled length of the
       * Newton step
       * @param dq_relative_error Relative error for different quotient
       * computation
       *
       * Line search parameters:
       *
       * @param maximum_beta_failures Maximum number of beta-condition failures
       *
       * Fixed point and Picard parameters:
       *
       * @param anderson_subspace_size Anderson acceleration subspace size
       */
      AdditionalData(const SolutionStrategy &strategy = linesearch,
                     const unsigned int maximum_non_linear_iterations = 200,
                     const double       function_tolerance            = 0.0,
                     const double       step_tolerance                = 0.0,
                     const bool         no_init_setup                 = false,
                     const unsigned int maximum_setup_calls           = 0,
                     const double       maximum_newton_step           = 0.0,
                     const double       dq_relative_error             = 0.0,
                     const unsigned int maximum_beta_failures         = 0,
                     const unsigned int anderson_subspace_size        = 0);

      /**
       * Add all AdditionalData() parameters to the given ParameterHandler
       * object. When the parameters are parsed from a file, the internal
       * parameters are automatically updated.
       *
       * The following parameters are declared:
       *
       * @code
       * set Function norm stopping tolerance       = 0
       * set Maximum number of nonlinear iterations = 200
       * set Scaled step stopping tolerance         = 0
       * set Solution strategy                      = linesearch
       * subsection Fixed point and Picard parameters
       *   set Anderson acceleration subspace size = 5
       * end
       * subsection Linesearch parameters
       *   set Maximum number of beta-condition failures = 0
       * end
       * subsection Newton parameters
       *   set Maximum allowable scaled length of the Newton step = 0
       *   set Maximum iterations without matrix setup            = 0
       *   set No initial matrix setup                            = false
       *   set Relative error for different quotient computation  = 0
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
       * Make sure that this class lives longer than `prm`. Undefined behavior
       * will occur if you destroy this class, and then parse a parameter file
       * using `prm`.
       */
      void
      add_parameters(ParameterHandler &prm);

      /**
       * The solution strategy to use. If you choose SolutionStrategy::newton
       * or SolutionStrategy::linesearch, you have to provide also the function
       * residual(). If you choose SolutionStrategy::picard or
       * SolutionStrategy::fixed_point, you have to provide also the function
       * iteration_function().
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
       * If set to zero, default values provided by KINSOL will be used.
       */
      double function_tolerance;

      /**
       * A scalar used as a stopping tolerance on the minimum
       * scaled step length.
       *
       * If set to zero, default values provided by KINSOL will be used.
       */
      double step_tolerance;

      /**
       * Whether an initial call to the preconditioner or Jacobian
       * setup function should be made or not.
       *
       * A call to this function is useful when solving a sequence of problems,
       * in which the final preconditioner or Jacobian value from one problem
       * is to be used initially for the next problem.
       */
      bool no_init_setup;

      /**
       * The maximum number of nonlinear iterations that can be
       * performed between calls to the setup_jacobian() function.
       *
       * If set to zero, default values provided by KINSOL will be used,
       * and in practice this often means that KINSOL will re-use a
       * Jacobian matrix computed in one iteration for later iterations.
       */
      unsigned int maximum_setup_calls;

      /**
       * The maximum allowable scaled length of the Newton step.
       *
       * If set to zero, default values provided by KINSOL will be used.
       */
      double maximum_newton_step;

      /**
       * The relative error in computing $F(u)$, which is used in the
       * difference quotient approximation to the Jacobian matrix when the user
       * does not supply a solve_jacobian_system_matrix() function.
       *
       * If set to zero, default values provided by KINSOL will be used.
       */
      double dq_relative_error;

      /**
       * The maximum number of beta-condition failures in the
       * linesearch algorithm. Only used if
       * strategy==SolutionStrategy::linesearch.
       */
      unsigned int maximum_beta_failures;

      /**
       * The size of the subspace used with Anderson acceleration
       * in conjunction with Picard or fixed-point iteration.
       *
       * If you set this to 0, no acceleration is used.
       */
      unsigned int anderson_subspace_size;
    };

    /**
     * Constructor. It is possible to fine tune the SUNDIALS KINSOL solver by
     * passing an AdditionalData() object that sets all of the solver
     * parameters.
     *
     * @param data KINSOL configuration data
     * @param mpi_comm MPI communicator
     */
    KINSOL(const AdditionalData &data     = AdditionalData(),
           const MPI_Comm &      mpi_comm = MPI_COMM_WORLD);

    /**
     * Destructor.
     */
    ~KINSOL();

    /**
     * Solve the non linear system. Return the number of nonlinear steps taken
     * to converge. KINSOL uses the content of `initial_guess_and_solution` as
     * initial guess, and stores the final solution in the same vector.
     */
    unsigned int
    solve(VectorType &initial_guess_and_solution);

    /**
     * A function object that users need to supply and that is intended to
     * reinitize the given vector to its correct size, block structure (if
     * block vectors are used), and MPI communicator (if the vector is
     * distributed across multiple processors using MPI), along with any
     * other properties necessary.
     */
    std::function<void(VectorType &)> reinit_vector;

    /**
     * A function object that users should supply and that is intended to
     * compute the residual `dst = F(src)`. This function is only used if the
     * SolutionStrategy::newton or SolutionStrategy::linesearch strategies
     * were selected.
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (KINSOL will try to change its internal
     * parameters and attempt a new solution step)
     * - <0: Unrecoverable error the computation will be aborted and an
     * assertion will be thrown.
     */
    std::function<int(const VectorType &src, VectorType &dst)> residual;

    /**
     * A function object that users should supply and that is intended to
     * compute the iteration function $G(u)$ for the fixed point and Picard
     * iteration. This function is only used if the
     * SolutionStrategy::fixed_point or SolutionStrategy::picard strategies
     * were selected.
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (KINSOL will try to change its internal
     * parameters and attempt a new solution step)
     * - <0: Unrecoverable error; the computation will be aborted and an
     * assertion will be thrown.
     */
    std::function<int(const VectorType &src, VectorType &dst)>
      iteration_function;

    /**
     * A function object that users may supply and that is intended to
     * prepare the linear solver for subsequent calls to
     * solve_jacobian_system().
     *
     * The job of setup_jacobian() is to prepare the linear solver for
     * subsequent calls to solve_jacobian_system(), in the solution of linear
     * systems $Ax = b$. The exact nature of this system depends on the
     * SolutionStrategy that has been selected.
     *
     * In the cases strategy = SolutionStrategy::newton or
     * SolutionStrategy::linesearch, $A$ is the Jacobian $J = \partial
     * F/\partial u$. If strategy = SolutionStrategy::picard, $A$ is the
     * approximate Jacobian matrix $L$. If strategy =
     * SolutionStrategy::fixed_point, then linear systems do not arise, and this
     * function is never called.
     *
     * The setup_jacobian() function may call a user-supplied function, or a
     * function within the linear solver module, to compute Jacobian-related
     * data that is required by the linear solver. It may also preprocess that
     * data as needed for solve_jacobian_system(), which may involve calling a
     * generic function (such as for LU factorization) or, more generally,
     * build preconditioners from the assembled Jacobian. In any case, the
     * data so generated may then be used whenever a linear system is solved.
     *
     * The point of this function is that
     * setup_jacobian() function is not called at every Newton iteration,
     * but only as frequently as the solver determines that it is appropriate
     * to perform the setup task. In this way, Jacobian-related data generated
     * by setup_jacobian() is expected to be used over a number of Newton
     * iterations. KINSOL determines itself when it is beneficial to regenerate
     * the Jacobian and associated information (such as preconditioners
     * computed for the Jacobian), thereby saving the effort to regenerate
     * the Jacobian matrix and a preconditioner for it whenever possible.
     *
     * @param current_u Current value of $u$
     * @param current_f Current value of $F(u)$ or $G(u)$
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (KINSOL will try to change its internal
     * parameters and attempt a new solution step)
     * - <0: Unrecoverable error the computation will be aborted and an
     * assertion will be thrown.
     */
    std::function<int(const VectorType &current_u, const VectorType &current_f)>
      setup_jacobian;

    /**
     * @deprecated Versions of SUNDIALS after 4.0 no longer provide all
     *   of the information necessary for this callback (see below). Use the
     *   `solve_with_jacobian` callback described below.
     *
     * A function object that users may supply and that is intended to solve
     * a linear system with the Jacobian matrix. This function will be called by
     * KINSOL (possibly several times) after setup_jacobian() has been called at
     * least once. KINSOL tries to do its best to call setup_jacobian() the
     * minimum number of times. If convergence can be achieved without updating
     * the Jacobian, then KINSOL does not call setup_jacobian() again. If, on
     * the contrary, internal KINSOL convergence tests fail, then KINSOL calls
     * setup_jacobian() again with updated vectors and coefficients so that
     * successive calls to solve_jacobian_systems() lead to better convergence
     * in the Newton process.
     *
     * If you do not specify a `solve_jacobian_system` or `solve_with_jacobian`
     * function, then only a fixed point iteration strategy can be used. Notice
     * that this may not converge, or may converge very slowly.
     *
     * A call to this function should store in `dst` the result of $J^{-1}$
     * applied to `rhs`, i.e., `J*dst = rhs`. It is the user's responsibility
     * to set up proper solvers and preconditioners inside this function
     * (or in the `setup_jacobian` callback above).
     *
     *
     * Arguments to the function are:
     *
     * @param[in] ycur The current $y$ vector for the current KINSOL
     * internal step. In the documentation above, this $y$ vector is generally
     * denoted by $u$.
     * @param[in] fcur The current value of the implicit right-hand side at
     * `ycur`, $f_I (t_n, ypred)$.
     * @param[in] rhs The system right hand side to solve for
     * @param[out] dst The solution of $J^{-1} * src$
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (KINSOL will try to change its internal
     * parameters and attempt a new solution step)
     * - <0: Unrecoverable error the computation will be aborted and an
     * assertion will be thrown.
     *
     * @warning Starting with SUNDIALS 4.1, SUNDIALS no longer provides the
     *   `ycur` and `fcur` variables -- only `rhs` is provided and `dst`
     *   needs to be returned. The first two arguments will therefore be
     *   empty vectors in that case. In practice, that means that one
     *   can no longer compute a Jacobian matrix for the current iterate
     *   within this function. Rather, this has to happen inside the
     *   `setup_jacobian` function above that receives this information.
     *   If it is important that the Jacobian corresponds to the *current*
     *   iterate (rather than a re-used Jacobian matrix that had been
     *   computed in a previous iteration and that therefore corresponds
     *   to a *previous* iterate), then you will also have to set the
     *   AdditionalData::maximum_newton_step variable to one, indicating
     *   that the Jacobian should be re-computed in every iteration.
     */
    DEAL_II_DEPRECATED_EARLY
    std::function<int(const VectorType &ycur,
                      const VectorType &fcur,
                      const VectorType &rhs,
                      VectorType &      dst)>
      solve_jacobian_system;

    /**
     * A function object that users may supply and that is intended to solve
     * a linear system with the Jacobian matrix. This function will be called by
     * KINSOL (possibly several times) after setup_jacobian() has been called at
     * least once. KINSOL tries to do its best to call setup_jacobian() the
     * minimum number of times. If convergence can be achieved without updating
     * the Jacobian, then KINSOL does not call setup_jacobian() again. If, on
     * the contrary, internal KINSOL convergence tests fail, then KINSOL calls
     * setup_jacobian() again with updated vectors and coefficients so that
     * successive calls to solve_jacobian_systems() lead to better convergence
     * in the Newton process.
     *
     * If you do not specify a `solve_with_jacobian` function, then only a
     * fixed point iteration strategy can be used. Notice that this may not
     * converge, or may converge very slowly.
     *
     * A call to this function should store in `dst` the result of $J^{-1}$
     * applied to `rhs`, i.e., `J*dst = rhs`. It is the user's responsibility
     * to set up proper solvers and preconditioners inside this function
     * (or in the `setup_jacobian` callback above). The function attached
     * to this callback is also provided with a tolerance to the linear solver,
     * indicating that it is not necessary to solve the linear system with
     * the Jacobian matrix exactly, but only to a tolerance that KINSOL will
     * adapt over time.
     *
     * Arguments to the function are:
     *
     * @param[in] rhs The system right hand side to solve for.
     * @param[out] dst The solution of $J^{-1} * src$.
     * @param[in] tolerance The tolerance with which to solve the linear system
     *   of equations.
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (KINSOL will try to change its internal
     * parameters and attempt a new solution step)
     * - <0: Unrecoverable error the computation will be aborted and an
     * assertion will be thrown.
     */
    std::function<
      int(const VectorType &rhs, VectorType &dst, const double tolerance)>
      solve_with_jacobian;

    /**
     * A function object that users may supply and that is intended to return a
     * vector whose components are the weights used by KINSOL to compute the
     * vector norm of the solution. The implementation of this function is
     * optional, and it is used only if implemented.
     *
     * The intent for this scaling factor is for problems in which the different
     * components of a solution have vastly different numerical magnitudes --
     * typically because they have different physical units and represent
     * different things. For example, if one were to solve a nonlinear Stokes
     * problem, the solution vector has components that correspond to velocities
     * and other components that correspond to pressures. These have different
     * physical units and depending on which units one chooses, they may have
     * roughly comparable numerical sizes or maybe they don't. To give just one
     * example, in simulations of flow in the Earth's interior, one has
     * velocities on the order of maybe ten centimeters per year, and pressures
     * up to around 100 GPa. If one expresses this in SI units, this corresponds
     * to velocities of around $0.000,000,003=3 \times 10^{-9}$ m/s, and
     * pressures around $10^9 \text{kg}/\text{m}/\text{s}^2$, i.e., vastly
     * different. In such cases, computing the $l_2$ norm of a solution-type
     * vector (e.g., the difference between the previous and the current
     * solution) makes no sense because the norm will either be dominated by the
     * velocity components or the pressure components. The scaling vector this
     * function returns is intended to provide each component of the solution
     * with a scaling factor that is generally chosen as as the inverse of a
     * "typical velocity" or "typical pressure" so that upon multiplication of a
     * vector component by the corresponding scaling vector component, one
     * obtains a number that is of order of magnitude of one (i.e., a reasonably
     * small multiple of one times the typical velocity/pressure). The KINSOL
     * manual states this as follows: "The user should supply values $D_u$,
     * which are diagonal elements of the scaling matrix such that $D_u U$ has
     * all components roughly the same magnitude when $U$ is close to a
     * solution".
     *
     * If no function is provided to a KINSOL object, then this is interpreted
     * as implicitly saying that all of these scaling factors should be
     * considered as one.
     */
    std::function<VectorType &()> get_solution_scaling;

    /**
     * A function object that users may supply and that is intended to return a
     * vector whose components are the weights used by KINSOL to compute the
     * vector norm of the function evaluation away from the solution. The
     * implementation of this function is optional, and it is used only if
     * implemented.
     *
     * The point of this function and the scaling vector it returns is similar
     * to the one discussed above for `get_solution_scaling`, except that it is
     * for a vector that scales the components of the function $F(U)$, rather
     * than the components of $U$, when computing norms. As above, if no
     * function is provided, then this is equivalent to using a scaling vector
     * whose components are all equal to one.
     */
    std::function<VectorType &()> get_function_scaling;

    /**
     * Handle KINSOL exceptions.
     */
    DeclException1(ExcKINSOLError,
                   int,
                   << "One of the SUNDIALS KINSOL internal functions "
                   << "returned a negative error code: " << arg1
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
     * KINSOL configuration data.
     */
    AdditionalData data;

    /**
     * KINSOL memory object.
     */
    void *kinsol_mem;

    /**
     * KINSOL solution vector.
     */
    N_Vector solution;

    /**
     * KINSOL solution scale.
     */
    N_Vector u_scale;

    /**
     * KINSOL f scale.
     */
    N_Vector f_scale;

    /**
     * MPI communicator. SUNDIALS solver runs happily in parallel.
     */
    MPI_Comm communicator;

    /**
     * Memory pool of vectors.
     */
    GrowingVectorMemory<VectorType> mem;
  };

} // namespace SUNDIALS


DEAL_II_NAMESPACE_CLOSE

#endif

#endif
