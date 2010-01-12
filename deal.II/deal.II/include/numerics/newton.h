//---------------------------------------------------------------------------
//    $Id: timestepping.cc 822 2010-01-11 08:57:30Z kanschat $
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__newton_step_control_h
#define __deal2__newton_step_control_h

#include <base/smartpointer.h>
#include <lac/solver_control.h>

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
 * than #threshold, the event #bad_derivative is submitted to
 * #inverse_derivative. It is up to this object to implement
 * reassembling accordingly.
 *
 * @author Guido Kanschat, 2006, 2010
 */
  template <class VECTOR>
  class Newton : public Operator<VECTOR>
  {
    public:
				       /**
					* Constructor, receiving the
					* application computing the
					* residual and solving the
					* linear problem.
					*/
      Newton (Operator<VECTOR>& residual, Operator<VECTOR>& inverse_derivative);
      
				       /**
					* Declare the parameter
					* applicable to Newton's method.
					*/
      void declare_parameters (ParameterHandler& param);
      
				       /**
					* Read the parameters.
					*/
      void initialize (ParameterHandler& param);
      
      				       /**
					* The actual Newton iteration.
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
      SmartPointer<Operator<VECTOR> > residual;

				       /**
					* The operator applying the
					* inverse derivative to the residual.
					*/
      SmartPointer<Operator<VECTOR> > inverse_derivative;
      
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
