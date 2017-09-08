//-----------------------------------------------------------
//
//    Copyright (C) 2017 by the deal.II authors
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

#ifndef dealii_sundials_ida_h
#define dealii_sundials_ida_h

#include <deal.II/base/config.h>
#ifdef DEAL_II_WITH_SUNDIALS

#include <deal.II/base/logstream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_view.h>
#include <deal.II/lac/vector_memory.h>

#include <ida/ida.h>
#include <ida/ida_spils.h>
#include <ida/ida_spgmr.h>
#include <ida/ida_spbcgs.h>
#include <ida/ida_sptfqmr.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN

// Shorthand notation for IDA error codes.
#define AssertIDA(code) Assert(code >= 0, ExcIDAError(code))

namespace SUNDIALS
{

  /**
   * Interface to SUNDIALS IDA (Implicit Differential-Algebraic) solver.
   *
   * The class IDAInterface is a wrapper to the Implicit
   * Differential-Algebraic solver which is a general purpose solver for
   * systems of Differential-Algebraic Equations (DAEs).
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
   *  - output_step;
   *  - get_local_tolerances;
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
   * where $y,\dot y$ are vectors in $\R^n$, $t$ is often the time (but can
   * also be a parametric quantity), and
   * $F:\R\times\R^n\times\R^n\rightarrow\R^n$. Such problem is solved
   * using Newton iteration augmented with a line search global
   * strategy. The integration method used in IDA is the variable-order,
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
   * application of the BDF method to the DAE system results in a nonlinear algebraic
   * system to be solved at each time step:
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
   * and $\alpha = \alpha_{n,0}/h_n$. It is worth metioning that the
   * scalar $\alpha$ changes whenever the step size or method order
   * changes.
   *
   * To provide a simple example, consider the following harmonic oscillator problem:
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
   *   y_0' & -y_1      & = 0 \\
   *   y_1' & + k^2 y_0 & = 0
   * \end{matrix}
   * \f]
   *
   * That is $F(y', y, t) = y' + A y = 0 $
   *
   * where
   * \f[
   * \begin{matrix}
   * 0 & -1 \\
   * k^2 &0
   * \end{matrix}
   * \f]
   *
   * and $y(0)=(0, k)$, $y'(0) = (k, 0)$
   *
   * The exact solution is $y_0(t) = \sin(k t)$, $y_1(t) = y_0'(t) = k \cos(k t)$,
   * $y_1'(t) = -k^2 \sin(k t)$
   *
   * The Jacobian to assemble is the following:  $J = \alpha I + A$
   *
   * This is achieved by the following snippet of code:
   * @code
   * typedef Vector<double> VectorType;
   *
   * VectorType y(2);
   * VectorType y_dot(2);
   *
   * FullMatrix<double> A(2,2);
   * A(0,1) = -1.0;
   * A(1,0) = kappa*kappa;
   *
   * FullMatrix<double> J(2,2);
   * FullMatrix<double> Jinv(2,2);
   *
   * double kappa = 1.0;
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
  template<typename VectorType=Vector<double> >
  class IDA
  {
  public:

    /**
     * Constructor. It is possible to fine tune the SUNDIALS IDA solver by tweaking some of
     * the input parameters.
     *
     * IDA is a Differential Algebraic solver. As such, it requires initial conditions also for
     * the first order derivatives. If you do not provide consistent initial conditions, (i.e.,
     * conditions for which F(y_dot(0), y(0), 0) = 0), you can ask SUNDIALS to compute initial
     * conditions for you by using the `ic_type` parameter at construction time.
     *
     * You have three options
     * -  none: do not try to make initial conditions consistent.
     * -  use_y_diff: compute the algebraic components of y and differential
     *    components of y_dot, given the differential components of y.
     *    This option requires that the user specifies differential and
     *    algebraic components in the function get_differential_components.
     * -  use_y_dot: compute all components of y, given y_dot.
     *
     * By default, this class assumes that all components are differential, and that you want
     * to solve a standard ode. In this case, the initial component type is set to `use_y_diff`,
     * so that the `y_dot` at time t=`initial_time` is computed by solving the nonlinear problem
     * $F(y_dot, y(t0), t0) = 0$ in the variable `y_dot`.
     *
     * Notice that a Newton solver is used for this computation. The Newton solver parameters
     * can be tweaked by acting on `ic_alpha` and `ic_max_iter`.
     *
     * If you reset the solver at some point, you may want to select a different computation
     * for the initial conditions after reset. Say, for example, that you have refined a grid,
     * and after transferring the solution to the new grid, the initial conditions are no longer
     * consistent. Then you can choose how these are made consistent, using the same three
     * options that you used for the initial conditions in `reset_type`.
     *
     * @param mpi_comm MPI communicator
     * @param initial_time Initial time
     * @param final_time Final time
     * @param initial_step_size Initial step size
     * @param min_step_size Minimum step size
     * @param abs_tol Absolute error tolerance
     * @param rel_tol Relative error tolerance
     * @param max_order Maximum BDF order
     * @param output_period Time interval between each output
     * @param ignore_algebraic_terms_for_errors Ignore algebraic terms for error computations
     * @param ic_type Initial condition type
     * @param reset_type Initial condition type after restart
     * @param ic_alpha Initial condition Newton parameter
     * @param ic_max_iter Initial condition Newton max iterations
     * @param max_non_linear_iterations Maximum number of nonlinear iterations
     * @param verbose Show output of time steps
     * @param use_local_tolerances Use local tolerances
     */
    IDA(const MPI_Comm mpi_comm = MPI_COMM_WORLD,
        const double &initial_time = 0.0,
        const double &final_time = 1.0,
        const double &initial_step_size = 1e-4,
        const double &min_step_size = 1e-6,
        const double &abs_tol = 1e-6,
        const double &rel_tol = 1e-5,
        const unsigned int &max_order = 5,
        const double &output_period = .1,
        const bool &ignore_algebraic_terms_for_errors = true,
        const std::string &ic_type = "none",
        const std::string &reset_type = "none",
        const double &ic_alpha = .33,
        const unsigned int &ic_max_iter = 5,
        const unsigned int &max_non_linear_iterations = 10,
        const bool &verbose = true,
        const bool &use_local_tolerances = false);

    /**
     * House cleaning.
     */
    ~IDA();

    /**
     * Add all internal parameters to the given ParameterHandler object. When
     * the paramaters are parsed from a file, the internal parameters are automatically
     * updated.
     *
     * The following parameters are declared:
     *
     * @code
     * set Absolute error tolerance                      = 0.000001
     * set Final time                                    = 1.000000
     * set Ignore algebraic terms for error computations = true
     * set Initial condition Newton max iterations       = 5
     * set Initial condition Newton parameter            = 0.330000
     * set Initial condition type                        = use_y_diff
     * set Initial condition type after restart          = use_y_dot
     * set Initial step size                             = 0.000100
     * set Initial time                                  = 0.000000
     * set Maximum number of nonlinear iterations        = 10
     * set Maximum order of BDF                          = 5
     * set Min step size                                 = 0.000001
     * set Relative error tolerance                      = 0.000010
     * set Show output of time steps                     = true
     * set Time units between each output                = 0.05
     * set Use local tolerances                          = false
     * @endcode
     *
     * These are one-to-one with the options you can pass at construction time.
     *
     * Those options are set as default value in the ParameterHandler object
     * `prm`.
     */
    virtual void add_parameters(ParameterHandler &prm);

    /**
     * Evolve. This function returns the final number of steps.
     */
    unsigned int solve_dae(VectorType &solution,
                           VectorType &solution_dot);

    /**
     * Clear internal memory and start with clean objects. This function is
     * called when the simulation start and when the user returns true to a
     * call to solver_should_restart().
     *
     * By default solver_should_restart returns false. If the user needs to
     * implement, for example, local adaptivity in space, he or she may assign
     * a different function to solver_should_restart() that performs all mesh
     * changes, transfers the solution and the solution dot to the new mesh, and
     * returns true.
     *
     * During reset(), both y and yp are checked for consistency, and according to
     * what was specified as ic_type (if t==initial_time) or reset_type (if
     * t>initial_time), yp, y, or both are modified to obtain a consistent set
     * of initial data.
     *
     * @param[in] t  The new starting time
     * @param[in] h  The new (tentative) starting time step
     * @param[in,out] y   The new (tentative) initial solution
     * @param[in,out] yp  The new (tentative) initial solution_dot
     */
    void reset(const double &t,
               const double &h,
               VectorType &y,
               VectorType &yp);

    /**
     * Reinit vector to have the right size, MPI communicator, etc.
     */
    std::function<void(VectorType &)> reinit_vector;

    /**
     * Compute residual.
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (IDAReinit will be called if this happens, and
     *       then last function will be attempted again
     * - <0: Unrecoverable error the computation will be aborted and an assertion
     *       will be thrown.
     */
    std::function<int(const double t,
                      const VectorType &y,
                      const VectorType &y_dot,
                      VectorType &res)> residual;

    /**
     * Compute Jacobian.
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (IDAReinit will be called if this happens, and
     *       then last function will be attempted again
     * - <0: Unrecoverable error the computation will be aborted and an assertion
     *       will be thrown.
     */
    std::function<int(const double t,
                      const VectorType &y,
                      const VectorType &y_dot,
                      const double alpha)> setup_jacobian;

    /**
     * Solve linear system.
     *
     * This function should return:
     * - 0: Success
     * - >0: Recoverable error (IDAReinit will be called if this happens, and
     *       then last function will be attempted again
     * - <0: Unrecoverable error the computation will be aborted and an assertion
     *       will be thrown.
     */
    std::function<int(const VectorType &rhs, VectorType &dst)> solve_jacobian_system;

    /**
     * Store solutions to file.
     */
    std::function<void (const double t,
                        const VectorType &sol,
                        const VectorType &sol_dot,
                        const unsigned int step_number)> output_step;

    /**
     * Evaluate wether the solver should be restarted (for example because the
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
    std::function<bool (const double t,
                        VectorType &sol,
                        VectorType &sol_dot)> solver_should_restart;

    /**
     * Return an index set containing the differential components.
     * Implementation of this function is optional. The default is to return a
     * complete index set. If your equation is also algebraic (i.e., it
     * contains algebraic constraings, or lagrange multipliers), you should
     * overwrite this function in order to return only the differential
     * components of your system.
     */
    std::function<IndexSet()> differential_components;

    /**
     * Return a vector whose components are the weights used by IDA to compute
     * the vector norm. The implementation of this function is optional, and it
     * is used only when `use_local_tolerances` is set to true at construction
     * time, or through the parameter file.
     *
     * If you do not overwrite this function and set `use_local_tolerances` to
     * true, an exception will be thrown when trying to start the time stepper.
     */
    std::function<VectorType&()> get_local_tolerances;

    /**
     * Set initial time equal to @p t disregarding what is written
     * in the parameter file.
     */
    void set_initial_time(const double &t);

    /**
     * Handle IDA exceptions.
     */
    DeclException1(ExcIDAError, int, << "One of the SUNDIALS IDA internal functions "
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
    double min_step_size;

    /**
     * Absolute error tolerance for adaptive time stepping.
     */
    double abs_tol;

    /**
     * Relative error tolerance for adaptive time stepping.
     */
    double rel_tol;

    /**
     * Maximum order of BDF.
     */
    unsigned int max_order;

    /**
     * Time period between each output.
     */
    double output_period;

    /**
     * Ignore algebraic terms for errors.
     */
    bool ignore_algebraic_terms_for_errors;

    /**
     * Type of initial conditions.
     *
     * IDA is a Differential Algebraic solver. As such, it requires initial conditions also for
     * the first order derivatives. If you do not provide consistent initial conditions, (i.e.,
     * conditions for which $F(y_dot(0), y(0), 0) = 0$), you can ask SUNDIALS to compute initial
     * conditions for you by using the `ic_type` parameter at construction time.
     *
     * You have three options
     * -  none: do not try to make initial conditions consistent.
     * -  use_y_diff: compute the algebraic components of y and differential
     *    components of y_dot, given the differential components of y.
     *    This option requires that the user specifies differential and
     *    algebraic components in the function get_differential_components.
     * -  use_y_dot: compute all components of y, given y_dot.
     *
     * Notice that you could in principle use this capabilities to solve for
     * stady state problems by setting y_dot to zero, and asking to compute
     * $y(0)$ that satisfies $F(0, y(0), 0) = 0$, however the nonlinear solver
     * used inside IDA may not be robust enough for complex problems with
     * several millions unknowns.
     */
    std::string ic_type;

    /**
     * Type of conditions to be used after a solver restart.
     *
     * If you do not have consistent initial conditions after a restart, (i.e.,
     * conditions for which F(y_dot(t_restart), y(t_restart), t_restart) = 0),
     * you can ask SUNDIALS to compute the new initial conditions for you by
     * using the `reset_type` parameter at construction time.
     *
     * You have three options
     * -  none: do not try to make initial conditions consistent.
     * -  use_y_diff: compute the algebraic components of y and differential
     *    components of y_dot, given the differential components of y.
     *    This option requires that the user specifies differential and
     *    algebraic components in the function get_differential_components.
     * -  use_y_dot: compute all components of y, given y_dot.
     */
    std::string reset_type;

    /**
     * Alpha to use in Newton method for IC calculation.
     */
    double ic_alpha;

    /**
     * Maximum number of iterations for Newton method in IC calculation.
     */
    unsigned ic_max_iter;

    /**
     * Maximum number of iterations for Newton method during time advancement.
     */
    unsigned int max_non_linear_iterations;

    /**
     * Show the progress of the time stepper.
     */
    bool verbose;

    /**
     * Use local tolerances when computing absolute tolerance.
     */
    bool use_local_tolerances;

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

#ifdef DEAL_II_WITH_MPI
    /**
     * MPI communicator. SUNDIALS solver runs happily in parallel.
     */
    MPI_Comm communicator;
#endif

    /**
     * Output stream. If run in parallel, only processor zero will
     * output verbose information about time steps.
     */
    ConditionalOStream pcout;

    /**
     * Memory pool of vectors.
     */
    GrowingVectorMemory<VectorType> mem;
  };

}


DEAL_II_NAMESPACE_CLOSE
#endif


#endif
