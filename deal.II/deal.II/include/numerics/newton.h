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

#ifndef __deal2__newton_h
#define __deal2__newton_h

#include <base/smartpointer.h>
#include <lac/solver_control.h>
#include <numerics/operator.h>

DEAL_II_NAMESPACE_OPEN

class ParameterHandler;

namespace Algorithms
{
/**
 * Operator class performing Newton's iteration with standard step
 * size control and adaptive matrix generation.
 *
 * This class performes a Newton iteration up to convergence
 * determined by #control. If after an update the norm of the residual
 * has become larger, then step size control is activated and the
 * update is subsequently divided by two until the residual actually
 * becomes smaller (or the minimal scaling factor determined by
 * #n_stepsize_iterations is reached).
 *
 * Since assembling matrices, depending on the implementation, tends
 * to be costly, this method applies an adaptive reassembling
 * strategy. Only if the reduction factor for the residual is more
 * than #threshold, the event Algorithms::bad_derivative is submitted to
 * #inverse_derivative. It is up to this object to implement
 * reassembling accordingly.
 *
 * <h3>Contents of the NamedData objects</h3>
 *
 * The only value used by the Newton method is the first vector in the
 * parameter <tt>out</tt> of operator()(). It serves as the start
 * vector of Newton's method and in the end contains the solution. All
 * other vectors of <tt>out</tt> are ignored by Newton's method and
 * its inner Operator objects. All vectors of <tt>in</tt> are forwarded to
 * the inner Operator objects, with additional information added as follows.
 *
 * When calling (*#residual)(), the NamedData <tt>in</tt> given to the
 * Newton iteration is prepended by a vector <tt>"Newton iterate"</tt>,
 * the current value of the Newton iterate, which can be used to
 * evaluate the residual at this point.
 *
 * For the call to (*#inverse_derivative), the vector <tt>"Newton
 * residual"</tt> is inserted before <tt>"Newton iterate"</tt>.
 *
 * @author Guido Kanschat, 2006, 2010
 */
  template <class VECTOR>
  class Newton : public Operator<VECTOR>
  {
    public:
				       /**
					* Constructor, receiving the
					* applications computing the
					* residual and solving the
					* linear problem, respectively.
					*/
      Newton (Operator<VECTOR>& residual, Operator<VECTOR>& inverse_derivative);
      
				       /**
					* Declare the parameters
					* applicable to Newton's method.
					*/
      void declare_parameters (ParameterHandler& param);
      
				       /**
					* Read the parameters.
					*/
      void initialize (ParameterHandler& param);

				       /**
					* Initialize the pointer 
                                        * data_out for debugging.
					*/
      void initialize (OutputOperator<VECTOR>& output);
      
      				       /**
					* The actual Newton
					* iteration. The initial value
					* is in <tt>out(0)</tt>, which
					* also contains the result
					* after convergence. Values in
					* <tt>in</tt> are not used by
					* Newton, but will be handed
					* down to the objects
					* #residual and #inverse_derivative.
					*/
      virtual void operator() (NamedData<VECTOR*>& out, const NamedData<VECTOR*>& in);
      
      virtual void notify(const Event&);

				       /**
					* Set the maximal residual
					* reduction allowed without
					* triggering assembling in the
					* next step. Return the
					* previous value.
					*/
      double threshold(double new_value);

				       /**
					* Control object for the
					* Newton iteration.
					*/
      ReductionControl control;
    private:
				       /**
					* The operator computing the residual.
					*/
      SmartPointer<Operator<VECTOR>, Newton<VECTOR> > residual;

				       /**
					* The operator applying the
					* inverse derivative to the residual.
					*/
      SmartPointer<Operator<VECTOR>, Newton<VECTOR> > inverse_derivative;

				       /**
					* The operator handling the output 
                                        * in case the debug_vectors is true.
                                        * Call the initialize function first.
					*/
      SmartPointer<OutputOperator<VECTOR>, Newton<VECTOR> > data_out;
      
				       /**
					* This flag is set by the
					* function assemble(),
					* indicating that the matrix
					* must be assembled anew upon
					* start.
					*/
      bool assemble_now;

				       /**
					* A flag used to decide how
					* many stepsize iteration
					* should be made. Default is
					* the original value of 21.
					*
					* Enter zero here to turn of
					* stepsize control.
					*
					* @note Controlled by
					* <tt>Stepsize iterations</tt> in
					* parameter file
					*/
      unsigned int n_stepsize_iterations;
      
				       /**
					* Threshold for re-assembling matrix.
					*
					* If the quotient of two
					* consecutive residuals is
					* smaller than this threshold,
					* the system matrix is not
					* assembled in this step.
					*
					* @note This parameter should be
					* adjusted to the residual gain
					* of the inner solver.
					*
					* The default values is zero,
					* resulting in reassembling in
					* every Newton step.
					*/
      double assemble_threshold;

    public:
				       /**
					* Print residual, update and
					* updated solution after each
					* step into file
					* <tt>Newton_NNN</tt>?
					*/
      bool debug_vectors;
				       /**
					* Write debug output to
					* #deallog; the higher the
					* number, the more output.
					*/
      unsigned int debug;	
  };
}

DEAL_II_NAMESPACE_CLOSE

#endif
