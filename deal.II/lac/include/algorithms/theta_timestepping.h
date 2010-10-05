//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__theta_timestepping_h
#define __deal2__theta_timestepping_h

#include <base/smartpointer.h>
#include <algorithms/operator.h>
#include <algorithms/timestep_control.h>

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
 * Assume that we want to solve the equation <i>u' + Au = 0</i> with a
 * step size <i>k</i>.  A step of the theta scheme can be written as
 *
 * @f[
 *   (M + \theta k A) u_{n+1} = (M - (1-\theta)k A) u_n.
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
 * autonomous).
 *
 * <h3>Usage of vectors in NamedData</h3>
 *
 * ThetaTimestepping uses NamedData for communicating vectors. With
 * outer or inner Operator objects. It does not use itself the input
 * vectors provided, but forwards them to the explicit and implicit
 * operators.
 *
 * The explicit Operator #op_explicit receives in its input in first
 * place the vector <tt>"Previous iterate"</tt>, which is the solution
 * value after the previous timestep. It is followed by all vectors
 * provided to ThetaTimestepping::operator() as input
 * argument. #op_explicit is supposed to write its result into the
 * first position of its output argument, labeled <tt>"Result"</tt>.
 *
 * The implicit Operator #op_implicit receives the result of
 * #op_explicit in its first input vector labeled <tt>"Previous
 * time"</tt>. It is followed by all vectors provided to
 * ThetaTimestepping::operator() as input argument. The output of
 * #op_implicit is directly written into the output argument given to
 * ThetaTimestepping.
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
 * us look at how they get used. First, let us define a matrix to be
 * used for our system and also an OutputOperator in order to write
 * the data of each timestep to a file.
 *
 * @skipline main
 * @until out.initialize
 *
 * Now we create objects for the implicit and explicit parts of the
 * steps as well as the ThetaTimestepping itself. Notice how the
 * TimestepData of ThetaTimestepping gets forwarded to the inner
 * operators. There are two different data objects, because the
 * timestep size is modified by #theta.
 *
 * @until set_output
 *
 * The next step is providing the vectors to be used. <tt>value</tt>
 * is filled with the initial value and is also the vector where the
 * solution at each timestep will be. Because the interface of
 * Operator has to be able to handle several vectors, we need to store
 * it in a NamedData object. Notice, that we need to create the
 * intermediate pointer <tt>p</tt>. If we would use
 * <code>&value</code> directly in the <code>add</code> function, the
 * resulting object would be constant.
 *
 * @until add
 *
 * Finally, we are ready to tell the solver, that we are looknig at
 * the initial timestep and run it.
 *
 * @until outdata
 * @skip Explicit::initialize
 *
 * Now we need to study the application of the implicit and explicit
 * operator. We assume that the pointer <code>matrix</code> points to
 * the matrix created in the main program, and that
 * <code>timestep_data</code> points to the correct data object of
 * ThetaTimestepping.
 *
 * @skipline void
 * @until vmult
 * @until }
 *
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
      ThetaTimestepping (Operator<VECTOR>& op_explicit,
			 Operator<VECTOR>& op_implicit);

      				       /**
					* The timestepping scheme. <tt>in</tt>
					* should contain the initial
					* value in first position. <tt>out</tt>
					*/
      virtual void operator() (NamedData<VECTOR*>& out, const NamedData<VECTOR*>& in);

				       /**
					* Register an event triggered
					* by an outer iteration.
					*/
      virtual void notify(const Event&);

				       /**
					* Define an operator which
					* will output the result in
					* each step. Note that no
					* output will be generated
					* without this.
					*/
      void set_output(OutputOperator<VECTOR>& output);
      
      void declare_parameters (ParameterHandler& param);
      void initialize (ParameterHandler& param);
				       /**
					* The current time in the
					* timestepping scheme.
					*/
      const double& current_time() const;
				       /**
					* The current step size.
					*/
      const double& step_size() const;
				       /**
					* The weight between implicit
					* and explicit part.
					*/
      const double& theta() const;

				       /**
					* The data handed to the
					* #op_explicit time stepping
					* operator.
					*
					* The time in here is the time
					* at the beginning of the
					* current step, the time step
					* is (1-#theta) times the
					* actual time step.
					*/
      const TimestepData& explicit_data() const;
      
				       /**
					* The data handed to the
					* #op_implicit time stepping
					* operator.
					*
					* The time in here is the time
					* at the beginning of the
					* current step, the time step
					* is #theta times the
					* actual time step.
					*/
      const TimestepData& implicit_data() const;

				       /**
					* Allow access to the control object.
					*/
      TimestepControl& timestep_control();
      
    private:
				       /**
					* The object controlling the
					* time step size and computing
					* the new time in each step.
					*/
      TimestepControl control;
      
				       /**
					* The control parameter theta in the
					* range <tt>[0,1]</tt>.
					*/
      double vtheta;
				       /**
					* Use adaptive #theta if
					* <tt>true</tt>.
					*/
      bool adaptive;

				       /**
					* The data for the explicit
					* part of the scheme.
					*/
      TimestepData d_explicit;
      
				       /**
					* The data for the implicit
					* part of the scheme.
					*/
      TimestepData d_implicit;
      
      
				       /**
					* The operator computing the
					* explicit part of the
					* scheme. This will receive in
					* its input data the value at
					* the current time with name
					* "Current time solution". It
					* should obtain the current
					* time and time step size from
					* explicit_data().
					*
					* Its return value is
					* <i>Mu+cAu</i>, where
					* <i>u</i> is the current
					* state vector, <i>M</i> the
					* mass matrix, <i>A</i> the
					* operator in space and
					* <i>c</i> is the time step
					* size in explicit_data().
					*/
      SmartPointer<Operator<VECTOR>, ThetaTimestepping<VECTOR> > op_explicit;
      
				       /**
					* The operator solving the
					* implicit part of the
					* scheme. It will receive in
					* its input data the vector
					* "Previous time". Information on the
					* timestep should be obtained
					* from implicit_data().
					*
					* Its return value is the
					* solution <i>u</i> of
					* <i>Mu-cAu=f</i>, where
					* <i>f</i> is the dual space
					* vector found in the
					* "Previous time" entry of the
					* input data, <i>M</i> the
					* mass matrix, <i>A</i> the
					* operator in space and
					* <i>c</i> is the time step
					* size in explicit_data().
					*/
      SmartPointer<Operator<VECTOR>, ThetaTimestepping<VECTOR> > op_implicit;

				       /**
					* The operator writing the
					* output in each time step
					*/
      SmartPointer<OutputOperator<VECTOR>, ThetaTimestepping<VECTOR> > output;
  };


  template <class VECTOR>
  inline
  const TimestepData&
  ThetaTimestepping<VECTOR>::explicit_data () const
  {
    return d_explicit;
  }
  

  template <class VECTOR>
  inline
  const TimestepData&
  ThetaTimestepping<VECTOR>::implicit_data () const
  {
    return d_implicit;
  }
  
  template <class VECTOR>
  inline
  void ThetaTimestepping<VECTOR>::set_output (OutputOperator<VECTOR>& out) 
  {
    output = &out;
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
