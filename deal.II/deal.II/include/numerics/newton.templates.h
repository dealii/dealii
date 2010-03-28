//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2006, 2007, 2008, 2009, 2010 by Guido Kanschat
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <numerics/newton.h>

#include <base/parameter_handler.h>
#include <base/logstream.h>
#include <lac/vector_memory.h>

#include <iomanip>


DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  template <class VECTOR>
  Newton<VECTOR>::Newton(Operator<VECTOR>& residual, Operator<VECTOR>& inverse_derivative)
		  :
		  residual(&residual), inverse_derivative(&inverse_derivative),
		  assemble_now(false),
		  n_stepsize_iterations(21),
		  assemble_threshold(0.),
		  debug_vectors(false),
		  debug(0)
  {}


  template <class VECTOR>
  void
  Newton<VECTOR>::declare_parameters(ParameterHandler& param)
  {
    param.enter_subsection("Newton");
    ReductionControl::declare_parameters (param);
    param.declare_entry("Assemble threshold", "0.", Patterns::Double());
    param.declare_entry("Stepsize iterations", "21", Patterns::Integer());
    param.declare_entry("Debug level", "0", Patterns::Integer());
    param.declare_entry("Debug vectors", "false", Patterns::Bool());
    param.leave_subsection();
  }

  template <class VECTOR>
  void
  Newton<VECTOR>::initialize (ParameterHandler& param)
  {
    param.enter_subsection("Newton");
    control.parse_parameters (param);
    assemble_threshold = param.get_double("Assemble threshold");
    n_stepsize_iterations = param.get_integer("Stepsize iterations");
    debug_vectors = param.get_bool("Debug vectors");
    param.leave_subsection ();
  }

  template <class VECTOR>
  void
  Newton<VECTOR>::initialize (OutputOperator<VECTOR>& output)
  {
    data_out = &output;  
  }

  template <class VECTOR>
  void
  Newton<VECTOR>::notify(const Event& e)
  {
    residual->notify(e);
    inverse_derivative->notify(e);
  }


  template <class VECTOR>
  void
  Newton<VECTOR>::operator() (NamedData<VECTOR*>& out, const NamedData<VECTOR*>& in)
  {
    Assert (out.size() == 1, ExcNotImplemented());
    deallog.push ("Newton");

    VECTOR& u = *out(0);

    if (debug>2)
      deallog << "u: " << u.l2_norm() << std::endl;

    GrowingVectorMemory<VECTOR> mem;
    typename VectorMemory<VECTOR>::Pointer Du(mem);
    typename VectorMemory<VECTOR>::Pointer res(mem);

    res->reinit(u);
    NamedData<VECTOR*> src1;
    NamedData<VECTOR*> src2;
    VECTOR* p = &u;
    src1.add(p, "Newton iterate");
    src1.merge(in);
    p = res;
    src2.add(p, "Newton residual");
    src2.merge(src1);
    NamedData<VECTOR*> out1;
    out1.add(p, "Residual");
    p = Du;
    NamedData<VECTOR*> out2;
    out2.add(p, "Update");

    unsigned int step = 0;
				     // fill res with (f(u), v)
    (*residual)(out1, src1);
    double resnorm = res->l2_norm();
    double old_residual = resnorm / assemble_threshold + 1;

    while (control.check(step++, resnorm) == SolverControl::iterate)
      {
					 // assemble (Df(u), v)
	if (resnorm/old_residual >= assemble_threshold)
	  inverse_derivative->notify (Events::bad_derivative);

	Du->reinit(u);
	try
	  {
	    (*inverse_derivative)(out2, src2);
	  }
	catch (SolverControl::NoConvergence& e)
	  {
	    deallog << "Inner iteration failed after "
		    << e.last_step << " steps with residual "
		    << e.last_residual << std::endl;
	  }

	if (debug_vectors)
	  {
	    NamedData<VECTOR*> out;
            VECTOR* p = &u; out.add(p, "solution");
	    p = Du; out.add(p, "update");
	    p = res; out.add(p, "residual");
            *data_out << step;
	    *data_out << out;
	  }

	u.add(-1., *Du);
	old_residual = resnorm;
	(*residual)(out1, src1);
	resnorm = res->l2_norm();

					 // Step size control
	unsigned int step_size = 0;
	while (resnorm >= old_residual)
	  {
	    ++step_size;
	    if (step_size > n_stepsize_iterations)
	      {
		deallog << "No smaller stepsize allowed!";
		break;
	      }
	    if (control.log_history())
	      deallog << "Trying step size: 1/" << (1<<step_size)
		      << " since residual was " << resnorm << std::endl;
	    u.add(1./(1<<step_size), *Du);
	    (*residual)(out1, src1);
	    resnorm = res->l2_norm();
	  }
      }
    deallog.pop();

				     // in case of failure: throw
				     // exception
    if (control.last_check() != SolverControl::success)
      throw SolverControl::NoConvergence (control.last_step(),
					  control.last_value());
				     // otherwise exit as normal
  }
}


DEAL_II_NAMESPACE_CLOSE
