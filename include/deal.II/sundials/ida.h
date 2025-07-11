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


#ifndef dealii_sundials_ida_h
#define dealii_sundials_ida_h

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

#  ifdef DEAL_II_SUNDIALS_WITH_IDAS
#    include <idas/idas.h>
#  else
#    include <ida/ida.h>
#  endif

#  include <deal.II/sundials/sundials_types.h>
#  include <deal.II/sundials/sunlinsol_wrapper.h>

#  include <boost/signals2.hpp>

#  include <nvector/nvector_serial.h>
#  include <sundials/sundials_config.h>
#  include <sundials/sundials_math.h>

#  include <memory>


DEAL_II_NAMESPACE_OPEN


namespace SUNDIALS
{
  /**
   * Interface to SUNDIALS Implicit Differential-Algebraic (IDA) solver.
   *
   * The class IDA is a wrapper to SUNDIALS Implicit Differential-Algebraic
   * solver which is a general purpose solver for systems of
   * Differential-Algebraic Equations (DAEs). Another class that can solve
   * this set of equations is PETScWrappers::TimeStepper.
   *
   * The user has to provide the implementation of the following std::functions:
   *  - reinit_vector;
   *  - residual;
   *  - setup_jacobian;
   *  - solve_with_jacobian;
   *
   * Optionally, also the following functions could be provided. By default
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
   *   J=\dfrac{\partial G}{\partial y}
   *   = \dfrac{\partial F}{\partial y} +
   *     \alpha \dfrac{\partial F}{\partial \dot y}\, ,
   * \f]
   *
   * and $\alpha = \alpha_{n,0}/h_n$. It is worth mentioning that the
   * scalar $\alpha$ changes whenever the step size or method order
   * changes.
   *
   *
   * <h3> A simple example: an ordinary differential equation </h3>
   *
   * To provide a simple example, consider the following harmonic oscillator
   * problem:
   * \f[ \begin{split}
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
   * That is, $F(y', y, t) = y' + A y = 0 $
   * where
   * \f[
   * A =
   * \begin{pmatrix}
   * 0 & -1 \\
   * k^2 &0
   * \end{pmatrix}
   * \f]
   * and $y(0)=(0, k)$, $y'(0) = (k, 0)$.
   *
   * The exact solution is $y_0(t) = \sin(k t)$, $y_1(t) = y_0'(t)
   * = k \cos(k t)$, $y_1'(t) = -k^2 \sin(k t)$.
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
   *                             VectorType &res)
   * {
   *   res = y_dot;
   *   A.vmult_add(res, y);
   * };
   *
   * time_stepper.setup_jacobian = [&](const double ,
   *                                   const VectorType &, // y
   *                                   const VectorType &, // y_dot
   *                                   const double alpha)
   * {
   *   J = A;
   *
   *   J(0,0) += alpha;
   *   J(1,1) += alpha;
   *
   *   Jinv.invert(J);
   * };
   *
   * time_stepper.solve_with_jacobian_system = [&](const VectorType &src,
   *                                               VectorType       &dst,
   *                                               const double)  // tolerance
   * {
   *   Jinv.vmult(dst,src);
   * };
   *
   * y[1] = kappa;
   * y_dot[0] = kappa;
   * time_stepper.solve_dae(y,y_dot);
   * @endcode
   *
   *
   * <h3> A differential algebraic equation (DAE) example </h3>
   *
   * A more interesting example is a situation where the form $F(y', y, t) = 0$
   * provides something genuinely more flexible than a typical ordinary
   * differential equation. Specifically, consider the equation
   * @f{align*}{
   *   u'(t) &= av(t),
   *   \\
   *   0 &= v(t) - u(t).
   * @f}
   * One can combine the two variables into $y(t) = [u(t), v(t)]^T$.
   * Here, one of the two variables does not have a time derivative. In
   * applications, this is often the case when one variable evolves in
   * time (here, $u(t)$) on its own time scale, and the other one finds
   * its value as a function of the former on a much faster time scale.
   * In the current context, we could of course easily eliminate $v(t)$
   * using the second equation, and would then just be left with the
   * equation
   * @f[
   *   u'(t) = au(t)
   * @f]
   * which has solution $u(t) = u(0)e^{at}$. But this is, in general, not
   * easily possible if the two variables are related by differential
   * operators. In fact, this happens quite frequently in application. Take,
   * for example, the time-dependent Stokes equations:
   * @f{align*}{
   *   \frac{\partial \mathbf u(\mathbf x,t)}{\partial t}
   *   - \nu \Delta \mathbf u(\mathbf x,t) + \nabla p(\mathbf x,t)
   *   &= \mathbf f(\mathbf x,t),
   *   \\
   *   \nabla \cdot \mathbf u(\mathbf x,t) &= 0.
   * @f}
   * Here, the fluid velocity $\mathbf u(\mathbf x,t)$ evolves over time,
   * and the pressure is always in equilibrium with the flow because the Stokes
   * equations are derived under the assumption that the speed of sound (at
   * which pressure perturbations propagate) is much larger than the fluid
   * velocity. As a consequence, there is no time derivative on the pressure
   * available in the equation, but unlike the simple model problem above, the
   * pressure can not easily be eliminated from the system. Similar situations
   * happen in step-21, step-31, step-32, step-43, and others, where a subset of
   * variables is always in instantaneous equilibrium with another set of
   * variables that evolves on a slower time scale.
   *
   * Another case where we *could* eliminate a variable but do not want to
   * is where that additional variable is introduced in the first place to work
   * around some other problem. As an example, consider the time dependent
   * version of the biharmonic problem we consider in step-47 (as well as some
   * later ones). The equations we would then be interested in would read
   * @f{align*}{
   *   \frac{\partial u(\mathbf x,t)}{\partial t} + \Delta^2 u(\mathbf x,t) &=
   *   f(\mathbf x,t).
   * @f}
   * As discussed in step-47, the difficulty is the presence of the fourth
   * derivatives. One way in which one can address this is by introducing
   * an auxiliary variable $v=\Delta u$ which would render the problem into
   * the following one that only ever has second derivatives which we know
   * how to deal with:
   * @f{align*}{
   *   \frac{\partial u(\mathbf x,t)}{\partial t} + \Delta v(\mathbf x,t) &=
   *   f(\mathbf x,t),
   *   \\
   *   v(\mathbf x,t)-\Delta u(\mathbf x,t) &= 0.
   * @f}
   * Here, the introduction of the additional variable was voluntary, and
   * could be undone, but we don't want that of course. Rather, we end
   * up with a differential-algebraic equation because the equations do
   * not have a time derivative for $v$.
   *
   * Rather than show how to solve the trivial (linear) case above, let us
   * instead consider the situation where we introduce another variable $v$ that
   * is related to $u$ by the nonlinear relationship $v=u^p$, $p\ge 1$:
   * @f{align*}{
   *   u'(t) &= a v(t)^{1/p},
   *   \\
   *   0 &= v(t) - u(t)^p.
   * @f}
   * We will impose initial conditions as
   * @f{align*}{
   *   u(0) &= 1 \\
   *   v(0) &= 1.
   * @f}
   * The problem continues to have the solution $u(t)=e^{at}$ with the
   * auxiliary variable satisfying $v(t)=[e^{at}]^p$. One would implement
   * all of this using the following little program where you have to recall
   * that
   * @f[
   *   F = \begin{pmatrix}u' -a v^{1/p} \\ -u^p + v \end{pmatrix}
   * @f]
   * and that the Jacobian we need to provide is
   * @f[
   *   J(\alpha) =
   *   = \dfrac{\partial F}{\partial y} +
   *     \alpha \dfrac{\partial F}{\partial \dot y}
   *   = \begin{pmatrix} \alpha && -av^{1/p-1}/p \\ -pu^{p-1} & 1 \end{pmatrix}
   * @f]
   *
   * All of this can be implemented using the following code:
   * @code
   *   const double a = 1.0;
   *   const double p = 1.5;
   *
   *   using VectorType = Vector<double>;
   *
   *   VectorType         y(2);
   *   VectorType         y_dot(2);
   *   FullMatrix<double> J(2, 2);
   *   FullMatrix<double> A(2, 2);
   *   FullMatrix<double> Jinv(2, 2);
   *
   *   SUNDIALS::IDA<VectorVector> time_stepper;
   *
   *   time_stepper.reinit_vector = [&](VectorType &v) {
   *     v.reinit(2);
   *   };
   *
   *   time_stepper.residual = [&](const double      t,
   *                               const VectorType &y,
   *                               const VectorType &y_dot,
   *                               VectorType       &res) {
   *     //  F(Y', Y, t) = [x' -a y^{1/p} ; -x^p + y]
   *     res    = 0;
   *     res[0] = y_dot[0] - a * std::pow(y[1], 1./p);
   *     res[1] = -std::pow(y[0], p) + y[1];
   *   };
   *
   *   time_stepper.setup_jacobian = [&](const double,
   *                                     const VectorType &y,
   *                                     const VectorType &,
   *                                     const double alpha) {
   *     // J = [alpha -ay^{1/p-1}/p ; -px^{p-1} 1]
   *     J(0, 0) = alpha;
   *     J(0, 1) = -a*std::pow(y[1], 1./p-1)/p;
   *     J(1, 0) = -p*std::pow(y[0], p-1);
   *     J(1, 1) = 1;
   *
   *     Jinv.invert(J);
   *   };
   *
   *   time_stepper.solve_with_jacobian =
   *     [&](const VectorType &src, VectorType &dst, const double) {
   *       Jinv.vmult(dst, src);
   *     };
   *
   *   // Provide initial values:
   *   y[0] = y[1] = 1;
   *   // Also provide initial derivatives. Note that
   *   //    v'(0) = d/dt[u^p](0) = p[u'(0)]^{p-1} = p a^{p-1}
   *   y_dot[0] = a;
   *   y_dot[1] = p*std::pow(a, p-1);
   *   time_stepper.solve_dae(y, y_dot);
   * @endcode
   * Note that in this code, we not only provide initial conditions for
   * $u$ and $v$, but also for $u'$ and $v'$. We can do this here because
   * we know what the exact solution is.
   *
   *
   * <h3> DAEs with missing initial conditions </h3>
   *
   * Whereas in the previous section, we were able to provide not only
   * initial values in the form of a vector for $y(0)$, but also for
   * $y'(0)$, this is not a common situation. For example, for the Stokes
   * equations mentioned above,
   * @f{align*}{
   *   \frac{\partial \mathbf u(\mathbf x,t)}{\partial t}
   *   - \nu \Delta \mathbf u(\mathbf x,t) + \nabla p(\mathbf x,t)
   *   &= \mathbf f(\mathbf x,t),
   *   \\
   *   \nabla \cdot \mathbf u(\mathbf x,t) &= 0,
   * @f}
   * one generally might have an initial velocity field for
   * $\mathbf u(\mathbf x,0)$, but typically one does not have an initial
   * pressure field $p(\mathbf x,0)$ nor either of these variables' time
   * derivatives at $t=0$.
   *
   * Fortunately, they can typically be computed via the relationship
   * $F(t,y,\dot y) = 0$. To illustrate how this can is done, let us
   * re-use the nonlinear example from the previous section:
   * @f{align*}{
   *   u'(t) &= a v(t)^{1/p},
   *   \\
   *   0 &= v(t) - u(t)^p.
   * @f}
   * If we now impose initial conditions for both variables, for
   * example
   * @f{align*}{
   *   u(0) &= 1 \\
   *   v(0) &= 1,
   * @f}
   * then the only change necessary is to create the time stepper via
   * @code
   *   SUNDIALS::IDA<VectorType>::AdditionalData data;
   *   data.ic_type = SUNDIALS::IDA<VectorType>::AdditionalData::use_y_diff;
   *   SUNDIALS::IDA<Vector<double>> time_stepper(data);
   * @endcode
   * and then we can run the program with the following at the end:
   * @code
   *   // Provide correct initial conditions y(0), but incorrect initial
   *   // derivatives y'(0):
   *   y[0] = y[1] = 1;  // correct
   *   y_dot[0]    = 0;  // wrong
   *   y_dot[1]    = 0;  // wrong
   *   time_stepper.solve_dae(y, y_dot);
   * @endcode
   * Here, IDA first compute $\dot y(0)$ before starting the time stepping
   * process.
   *
   * In many applications, however, one does not even have a complete set of
   * initial conditions -- e.g., in the Stokes equations above, one generally
   * only has initial values for the velocity, but not the pressure. IDA can
   * also compute these, but for that it needs to know which components of the
   * solution vector have differential equations attached to them -- i.e., for
   * which components a time derivative appears in $F(t,y,\dot y)$. This is
   * not difficult to do -- we only have to add the following block where a
   * lambda function returns an IndexSet that describes which variables
   * are "differential" (included in the index set) and which are not (not
   * included in the index set):
   * @code
   *     time_stepper.differential_components = []() {
   *     IndexSet x(2);
   *     x.add_index(0);
   *     return x;
   *   };
   *
   *   y[0]     = 1;   // correct
   *   y[1]     = 42;  // wrong
   *   y_dot[0] = 0;   // wrong
   *   y_dot[1] = 0;   // wrong
   *   time_stepper.solve_dae(y, y_dot);
   * @endcode
   * With these modifications, IDA correctly computes the solutions $u(t)$
   * and $v(t)$.
   *
   * A word of caution, however: All of this solving for components of
   * $y(0)$ and $y'(t)$ costs time and accuracy. If you *can* provide initial
   * conditions, you should; if you can't, they have to be numerically
   * approximated and will be close but not exact. In the examples above,
   * if all initial conditions $y(0),\dot y(0)$ are provided, IDA computes
   * the solution $y(10)=e^{10}\approx 22,000$ to an absolute accuracy of
   * around $3\cdot 10^{-5}$ (i.e., to a relative tolerance of better than
   * $10^{-8}$). If you only provide $y(0)$ correctly, the absolute error
   * is about twice as large, around $6\cdot 10^{-5}$. If one also omits
   * providing the initial value for the second component of $y(0)$ (the
   * non-differential component $v(0)$), the error goes up to $5\cdot 10^{-4}$.
   * That's not bad, but the trend is clear. In practice, one can control
   * the accuracy of the required solves for initial conditions by
   * setting the appropriate flags in the AdditionalData object passed to
   * the constructor.
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
       * consistent initial conditions, (i.e., conditions for which $F(\dot
       * y(0), y(0), 0) = 0)$, you can ask SUNDIALS to compute initial
       * conditions for you by specifying InitialConditionCorrection for the
       * initial conditions both at the `initial_time` (`ic_type`) and after a
       * reset has occurred (`reset_type`).
       */
      enum InitialConditionCorrection
      {
        /**
         * Do not try to make initial conditions consistent.
         */
        none = 0,

        /**
         * Compute the algebraic components of $y$ and differential
         * components of $\dot y$, given the differential components of $y$.
         * This option requires that the user specifies differential and
         * algebraic components in the function
         * IDA::differential_components().
         */
        use_y_diff = 1,

        /**
         * Compute all components of $y$, given $\dot y$.
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
       * @param ls_norm_factor Converting factor from the integrator tolerance
       * to the linear solver tolerance
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
        const double       ls_norm_factor                = 0,
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
        , ls_norm_factor(ls_norm_factor)
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
       * Make sure that this class lives longer than `prm`. Undefined behavior
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
          "    algebraic components in the function differential_components().\n"
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
          "    algebraic components in the function differential_components().\n"
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
        prm.add_parameter(
          "Factor to use when converting from the integrator tolerance to the linear solver tolerance",
          ls_norm_factor);
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

      /**
       * Factor to use when converting from the integrator tolerance to the
       * linear solver tolerance
       */
      double ls_norm_factor;
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
     *    algebraic components in the function differential_components().
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
     * @param data IDA configuration data
     *
     * @note With SUNDIALS 6 and later this constructor sets up logging
     * objects to only work on the present processor (i.e., results are only
     * communicated over MPI_COMM_SELF).
     */
    IDA(const AdditionalData &data = AdditionalData());

    /**
     * Constructor.
     *
     * @param data IDA configuration data. Same as the other constructor.
     * @param mpi_comm MPI Communicator over which logging operations are
     * computed. Only used in SUNDIALS 6 and newer.
     */
    IDA(const AdditionalData &data, const MPI_Comm mpi_comm);

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
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, IDA can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
     */
    std::function<void(const double      t,
                       const VectorType &y,
                       const VectorType &y_dot,
                       VectorType       &res)>
      residual;

    /**
     * Compute Jacobian. This function is called by IDA any time a Jacobian
     * update is required. The user should compute the Jacobian (or update all
     * the variables that allow the application of the Jacobian). This function
     * is called by IDA once, before any call to solve_with_jacobian().
     *
     * The Jacobian $J$ should be a (possibly inexact) computation of
     * \f[
     *   J=\dfrac{\partial G}{\partial y} = \dfrac{\partial F}{\partial y} +
     *  \alpha \dfrac{\partial F}{\partial \dot y}.
     * \f]
     *
     * If the user uses a matrix based computation of the Jacobian, then this
     * is the right place where an assembly routine should be called to
     * assemble both a matrix and a preconditioner for the Jacobian system.
     * Subsequent calls (possibly more than one) to solve_with_jacobian() can
     * assume that this function has been called at least once.
     *
     * Notice that no assumption is made by this interface on what the user
     * should do in this function. IDA only assumes that after a call to
     * setup_jacobian() it is possible to call solve_with_jacobian() to obtain a
     * solution $x$ to the system $J x = b$.
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, IDA can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
     */
    std::function<void(const double      t,
                       const VectorType &y,
                       const VectorType &y_dot,
                       const double      alpha)>
      setup_jacobian;

    /**
     * Solve the Jacobian linear system up to a specified tolerance. This
     * function will be called by IDA (possibly several times) after
     * setup_jacobian() has been called at least once. IDA tries to do its best
     * to call setup_jacobian() the minimum number of times. If convergence can
     * be achieved without updating the Jacobian, then IDA does not call
     * setup_jacobian() again. If, on the contrary, internal IDA convergence
     * tests fail, then IDA calls again setup_jacobian() with updated vectors
     * and coefficients so that successive calls to
     * solve_with_jacobian() lead to better convergence in the
     * Newton process.
     *
     * The Jacobian $J$ should be (an approximation of) the system Jacobian
     * \f[
     *   J=\dfrac{\partial G}{\partial y} = \dfrac{\partial F}{\partial y} +
     *  \alpha \dfrac{\partial F}{\partial \dot y}.
     * \f]
     *
     * Arguments to the function are:
     *
     * @param[in] rhs The system right hand side to solve for.
     * @param[out] dst The solution of $J^{-1} * src$.
     * @param[in] tolerance The tolerance with which to solve the linear system
     *   of equations.
     *
     * A call to this function should store in `dst` the result of $J^{-1}$
     * applied to `src`, i.e., the solution of the linear system `J*dst = src`.
     * It is the user's responsibility to set up proper solvers and
     * preconditioners either inside this function, or already within the
     * `setup_jacobian()` function. (The latter is, for example, what the
     * step-77 program does: All expensive operations happen in
     * `setup_jacobian()`, given that that function is called far less often
     * than the current one.)
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, IDA can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
     */
    std::function<
      void(const VectorType &rhs, VectorType &dst, const double tolerance)>
      solve_with_jacobian;

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
                       const VectorType  &sol,
                       const VectorType  &sol_dot,
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
     *
     * @note This variable represents a
     * @ref GlossUserProvidedCallBack "user provided callback".
     * See there for a description of how to deal with errors and other
     * requirements and conventions. In particular, IDA can deal
     * with "recoverable" errors in some circumstances, so callbacks
     * can throw exceptions of type RecoverableUserCallbackError.
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
     *
     * When running in parallel, every process will call this function
     * independently, and synchronization will happen at the end of the
     * initialization setup to communicate what components are local. Make sure
     * you only return the locally owned (or locally relevant) components, in
     * order to minimize communication between processes.
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
                   << "One of SUNDIALS IDA's internal functions "
                   << "returned an error code: " << arg1
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
    const AdditionalData data;

    /**
     * IDA memory object.
     */
    void *ida_mem;

#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    /**
     * A context object associated with the IDA solver.
     */
    SUNContext ida_ctx;
#  endif

    /**
     * MPI communicator. Only used for SUNDIALS' logging routines - the actual
     * solve routines will use the communicator provided by the vector class.
     */
    MPI_Comm mpi_communicator;

    /**
     * Memory pool of vectors.
     */
    GrowingVectorMemory<VectorType> mem;

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
} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SUNDIALS

#endif
