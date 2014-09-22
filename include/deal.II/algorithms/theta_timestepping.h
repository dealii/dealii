// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#ifndef __deal2__theta_timestepping_h
#define __deal2__theta_timestepping_h

#include <deal.II/base/smartpointer.h>
#include <deal.II/algorithms/operator.h>
#include <deal.II/algorithms/timestep_control.h>

DEAL_II_NAMESPACE_OPEN

class ParameterHandler;

namespace Algorithms
{
  /**
   * A little structure, gathering the size of a timestep and the
   * current time. Time stepping schemes can use this to provide time
   * step information to the classes actually performing a single step.
   *
   * The definition of what is considered "current time" depends on the
   * scheme. For an explicit scheme, this is the time at the beginning
   * of the step. For an implicit scheme, it is usually the time at the
   * end.
   */
  struct TimestepData
  {
/// The current time
    double time;
/// The current step size times something
    double step;
  };

  /**
   * Application class performing the theta timestepping scheme.
   *
   * The theta scheme is an abstraction of implicit and explicit Euler
   * schemes, the Crank-Nicholson scheme and linear combinations of
   * those. The choice of the actual scheme is controlled by the
   * parameter #theta as follows.
   * <ul>
   * <li> #theta=0: explicit Euler scheme
   * <li> #theta=1: implicit Euler scheme
   * <li> #theta=½: Crank-Nicholson scheme
   * </ul>
   *
   * For fixed #theta, the Crank-Nicholson scheme is the only second
   * order scheme. Nevertheless, further stability may be achieved by
   * choosing #theta larger than ½, thereby introducing a first order
   * error term. In order to avoid a loss of convergence order, the
   * adaptive theta scheme can be used, where <i>#theta=½+c dt</i>.
   *
   * Assume that we want to solve the equation <i>u' + F(u) = 0</i> with a
   * step size <i>k</i>.  A step of the theta scheme can be written as
   *
   * @f[
   *   M u_{n+1} + \theta k F(u_{n+1})  = M u_n - (1-\theta)k F(u_n).
   * @f]
   *
   * Here, <i>M</i> is the mass matrix. We see, that the right hand side
   * amounts to an explicit Euler step with modified step size in weak
   * form (up to inversion of M). The left hand side corresponds to an
   * implicit Euler step with modified step size (right hand side
   * given). Thus, the implementation of the theta scheme will use two
   * Operator objects, one for the explicit, one for the implicit
   * part. Each of these will use its own TimestepData to account for
   * the modified step sizes (and different times if the problem is not
   * autonomous). Note that once the explicit part has been computed,
   * the left hand side actually constitutes a linear or nonlinear
   * system which has to be solved.
   *
   * <h3>Usage AnyData</h3>
   *
   * ThetaTimestepping uses AnyData for communicating vectors and time
   * step information. With
   * outer or inner Operator objects. It does not use itself the input
   * vectors provided, but forwards them to the explicit and implicit
   * operators.
   *
   * <h4>Vector data</h4>
   *
   * The explicit Operator #op_explicit receives in its input in first
   * place the vector "Previous iterate", which is the solution
   * value after the previous timestep. It is followed by all vectors
   * provided to ThetaTimestepping::operator() as input
   * argument. #op_explicit is supposed to write its result into the
   * first position of its output argument, labeled "Result".
   *
   * The implicit Operator #op_implicit receives the result of
   * #op_explicit in its first input vector labeled "Previous
   * time". It is followed by all vectors provided to
   * ThetaTimestepping::operator() as input argument. The output of
   * #op_implicit is directly written into the output argument given to
   * ThetaTimestepping.
   *
   * <h4>Scalar data</h4>
   *
   * Since the introduction of AnyData, ThetaTimestepping is able to
   * communicate the current time step information through AnyData as
   * well. Therefore, the AnyData objects handed as input to
   * #op_explicit and #op_implicit contain two entries of type
   * `const double*` named "Time" and "Timestep". Note
   * that "Time" refers to the time at the beginning of the current
   * step for #op_explicit and at the end for #op_implicit,
   * respectively.
   *
   * <h3>Usage of ThetaTimestepping</h3>
   *
   * The use ThetaTimestepping is more complicated than for instance
   * Newton, since the inner operators will usually need to access the
   * TimeStepData. Thus, we have a circular dependency of information,
   * and we include the following example for its use. It can be found
   * in <tt>examples/doxygen/theta_timestepping.cc</tt>
   *
   * @dontinclude theta_timestepping.cc
   *
   * First, we define the two operators used by ThetaTimestepping and
   * call them <code>Implicit</code> and <code>Explicit</code>. They
   * both share the public interface of Operator, and additionally
   * provide storage for the matrices to be used and a pointer to
   * TimestepData. Note that we do not use a SmartPointer here, since
   * the TimestepData will be destroyed before the operator.
   *
   * @skip class Explicit
   * @until End of declarations
   *
   * These operators will be implemented after the main program. But let
   * us look first at how they get used. First, let us define a matrix to be
   * used for our system and also an OutputOperator in order to write
   * the data of each timestep to a file.
   *
   * @skipline main
   * @until out.initialize
   *
   * Now we create objects for the implicit and explicit parts of the
   * steps as well as the ThetaTimestepping itself. We initialize the
   * timestepping with the output operator in order to be able to see
   * the output in every step.
   *
   * @until set_output
   *
   * The next step is providing the vectors to be used. <tt>value</tt>
   * is filled with the initial value and is also the vector where the
   * solution at each timestep will be. Because the interface of
   * Operator has to be able to handle several vectors, we need to store
   * it in an AnyData object. Since our problem has no additional
   * parameters, the input AnyData object remains empty.
   *
   * @until add
   *
   * Finally, we are ready to tell the solver, that we are starting at
   * the initial timestep and run it.
   *
   * @until }
   *
   * First the constructor, which simply copies the system matrix into
   * the member pointer for later use.
   *
   * @skip Explicit::
   * @until }
   *
   * Now we need to study the application of the implicit and explicit
   * operator. We assume that the pointer <code>matrix</code> points to
   * the matrix created in the main program (the constructor did this
   * for us). Here, we first get the time step size from the AnyData
   * object that was provided as input. Then, if we are in the first
   * step or if the timestep has changed, we fill the local matrix
   * $m$, such that with the given matrix $M$, it becomes
   * \f[
   * m = I - \Delta t M.
   * \f]
   * After we have worked off the notifications, we clear them, such
   * that the matrix is only generated when necessary.
   *
   * @skipline void
   * @until clear
   *
   * Now we multiply the input vector with the new matrix and store on output.
   *
   * @until }
   * The code for the implicit operator is almost the same, except
   * that we change the sign in front of the timestep and use the
   * inverse of t he matrix.
   *
   * @until vmult
   * @until }
   * @author Guido Kanschat
   * @date 2010
   */
  template <class VECTOR>
  class ThetaTimestepping : public Operator<VECTOR>
  {
  public:
    /**
     * Constructor, receiving the
     * two operators stored in
     * #op_explicit and #op_implicit. For
     * their meening, see the
     * description of those variables.
     */
    ThetaTimestepping (Operator<VECTOR> &op_explicit,
                       Operator<VECTOR> &op_implicit);

    /**
     * The timestepping scheme.
     *
     * @param in is ignored by
     * ThetaTimestepping, but is merged into the AnyData objects used
     * as input for the operators #op_explicit and #op_implicit.
     *
     * @param out in its first argument must contain a pointer to a
     * `VECTOR`, which contains the initial value when the operator is
     * called. It contains the final value when the operator returns.
     */
    virtual void operator() (AnyData &out, const AnyData &in);

    /**
     * @deprecated Use  the function with AnyData arguments
     */
    virtual void operator() (NamedData<VECTOR *> &out, const NamedData<VECTOR *> &in);

    /**
     * Register an event triggered
     * by an outer iteration.
     */
    virtual void notify(const Event &);

    /**
     * Define an operator which will output the result in each
     * step. Note that no output will be generated without this.
     */
    void set_output(OutputOperator<VECTOR> &output);

    /**
     * Declare parameters in a parameter handler.
     */
    static void declare_parameters (ParameterHandler &param);

    /**
     * Read the parameters in the ParameterHandler.
     */
    void parse_parameters (ParameterHandler &param);

    /**
     * @deprecated Use parse_parameters().
     */
    void initialize (ParameterHandler &param) DEAL_II_DEPRECATED;
    /**
     * The current time in the
     * timestepping scheme.
     */
    double current_time() const;
    /**
     * The current step size.
     */
    double step_size() const;
    /**
     * The weight between implicit and explicit part.
     */
    double theta() const;

    /**
     * Set a new weight and return the old
     */
    double theta(double new_theta);

    /**
     * The data handed to the #op_explicit time stepping operator.
     *
     * The time in here is the time at the beginning of the current
     * step, the time step is (1-#theta) times the actual time step.
     */
    const TimestepData &explicit_data() const;

    /**
     * The data handed to the #op_implicit time stepping operator.
     *
     * The time in here is the time at the beginning of the current
     * step, the time step is #theta times the actual time step.
     */
    const TimestepData &implicit_data() const;

    /**
     * Allow access to the control object.
     */
    TimestepControl &timestep_control();

  private:
    /**
     * The object controlling the time step size and computing the new
     * time in each step.
     */
    TimestepControl control;

    /**
     * The control parameter theta in the range <tt>[0,1]</tt>. It
     * defaults to 0.5.
     */
    double vtheta;
    /**
     * Use adaptive #theta if <tt>true</tt>. Not yet implemented.
     */
    bool adaptive;

    /**
     * The data for the explicit part of the scheme.
     */
    TimestepData d_explicit;

    /**
     * The data for the implicit part of the scheme.
     */
    TimestepData d_implicit;


    /**
     * The operator computing the explicit part of the scheme. This
     * will receive in its input data the value at the current time
     * with name "Current time solution". It should obtain the current
     * time and time step size from explicit_data().
     *
     * Its return value is $ Mu+cF(u) $, where $u$ is the
     * current state vector, $M$ the mass matrix, $F$ the
     * operator in space and $c$ is the adjusted
     * time step size $(1-\theta) \Delta t$.
     */
    SmartPointer<Operator<VECTOR>, ThetaTimestepping<VECTOR> > op_explicit;

    /**
     * The operator solving the implicit part of the scheme. It will
     * receive in its input data the vector "Previous
     * time". Information on the timestep should be obtained from
     * implicit_data().
     *
     * Its return value is the solution <i>u</i> of <i>Mu-cF(u)=f</i>,
     * where <i>f</i> is the dual space vector found in the "Previous
     * time" entry of the input data, <i>M</i> the mass matrix,
     * <i>F</i> the operator in space and <i>c</i> is the adjusted
     * time step size $ \theta \Delta t$
     */
    SmartPointer<Operator<VECTOR>, ThetaTimestepping<VECTOR> > op_implicit;

    /**
     * The operator writing the output in each time step
     */
    SmartPointer<OutputOperator<VECTOR>, ThetaTimestepping<VECTOR> > output;
  };


  template <class VECTOR>
  inline
  const TimestepData &
  ThetaTimestepping<VECTOR>::explicit_data () const
  {
    return d_explicit;
  }


  template <class VECTOR>
  inline
  const TimestepData &
  ThetaTimestepping<VECTOR>::implicit_data () const
  {
    return d_implicit;
  }


  template <class VECTOR>
  inline
  TimestepControl &
  ThetaTimestepping<VECTOR>::timestep_control ()
  {
    return control;
  }

  template <class VECTOR>
  inline
  void ThetaTimestepping<VECTOR>::set_output (OutputOperator<VECTOR> &out)
  {
    output = &out;
  }


  template <class VECTOR>
  inline
  double ThetaTimestepping<VECTOR>::theta () const
  {
    return vtheta;
  }


  template <class VECTOR>
  inline
  double ThetaTimestepping<VECTOR>::theta (double new_theta)
  {
    const double tmp = vtheta;
    vtheta = new_theta;
    return tmp;
  }


  template <class VECTOR>
  inline
  double ThetaTimestepping<VECTOR>::current_time () const
  {
    return control.now();
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
