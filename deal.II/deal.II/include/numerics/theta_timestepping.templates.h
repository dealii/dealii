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

#include <numerics/theta_timestepping.h>

#include <base/parameter_handler.h>
#include <lac/vector_memory.h>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  template <class VECTOR>
  ThetaTimestepping<VECTOR>::ThetaTimestepping (Operator<VECTOR>& e, Operator<VECTOR>& i)
		  : op_explicit(&e), op_implicit(&i)
  {}


  template <class VECTOR>
  void
  ThetaTimestepping<VECTOR>::notify(const Event& e)
  {
    op_explicit->notify(e);
    op_implicit->notify(e);
  }

  template <class VECTOR>
  void
  ThetaTimestepping<VECTOR>::declare_parameters(ParameterHandler& param)
  {
    param.enter_subsection("ThetaTimestepping");
    TimestepControl::declare_parameters (param);
    param.declare_entry("Theta", ".5", Patterns::Double());
    param.declare_entry("Adaptive", "false", Patterns::Bool());
    param.leave_subsection();
  }

  template <class VECTOR>
  void
  ThetaTimestepping<VECTOR>::initialize (ParameterHandler& param)
  {
    param.enter_subsection("ThetaTimestepping");
    control.parse_parameters (param);
    vtheta = param.get_double("Theta");
    adaptive = param.get_bool("Adaptive");
    param.leave_subsection ();
  }


  template <class VECTOR>
  void
  ThetaTimestepping<VECTOR>::operator() (NamedData<VECTOR*>& out, const NamedData<VECTOR*>& in)
  {
    Assert(!adaptive, ExcNotImplemented());

    deallog.push ("Theta");
    GrowingVectorMemory<VECTOR> mem;
    typename VectorMemory<VECTOR>::Pointer aux(mem);
    aux->reinit(*out(0));

    control.restart();

    d_explicit.time = control.now();

				     // The data used to compute the
				     // vector associated with the old
				     // timestep
    VECTOR* p = out(0);
    NamedData<VECTOR*> src1;
    src1.add(p, "Previous time");
    src1.merge(in);
    p = aux;
    NamedData<VECTOR*> out1;
    out1.add(p, "Result");
				     // The data provided to the inner
				     // solver
    NamedData<VECTOR*> src2;
    src2.add(p, "Previous time data");
    src2.merge(in);

    if (output != 0)
      (*output) << 0U << out;

    for (unsigned int count = 1; d_explicit.time < control.final(); ++count)
      {
	const bool step_change = control.advance();
	d_implicit.time = control.now();
	d_explicit.step = (1.-vtheta)*control.step();
	d_implicit.step = vtheta*control.step();
	deallog << "Time step:" << d_implicit.time << std::endl;

	op_explicit->notify(Events::new_time);
	op_implicit->notify(Events::new_time);
	if (step_change)
	  {
	    op_explicit->notify(Events::new_timestep_size);
	    op_implicit->notify(Events::new_timestep_size);
	  }

					 // Compute
					 // (I + (1-theta)dt A) u
	(*op_explicit)(out1, src1);
	(*op_implicit)(out, src2);

	if (output != 0 && control.print())
	  (*output) << count << out;

	d_explicit.time = control.now();
      }
    deallog.pop();
  }
}

DEAL_II_NAMESPACE_CLOSE
