//---------------------------------------------------------------------------
//    $Id: time_dependent.all_dimensions.cc 21627 2010-08-09 05:10:10Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2005, 2006, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <numerics/time_dependent.h>
#include <base/memory_consumption.h>
#include <base/thread_management.h>

#include <functional>
#include <algorithm>
#include <numeric>

DEAL_II_NAMESPACE_OPEN

TimeDependent::TimeSteppingData::TimeSteppingData (const unsigned int look_ahead,
						   const unsigned int look_back)
		:
		look_ahead (look_ahead),
		look_back (look_back)
{}


TimeDependent::TimeDependent (const TimeSteppingData &data_primal,
			      const TimeSteppingData &data_dual,
			      const TimeSteppingData &data_postprocess):
		sweep_no (numbers::invalid_unsigned_int),
		timestepping_data_primal (data_primal),
		timestepping_data_dual (data_dual),
		timestepping_data_postprocess (data_postprocess)
{}


TimeDependent::~TimeDependent ()
{
  while (timesteps.size() != 0)
    delete_timestep (0);
}


void
TimeDependent::insert_timestep (const TimeStepBase *position,
				TimeStepBase       *new_timestep)
{
  Assert ((std::find(timesteps.begin(), timesteps.end(), position) != timesteps.end()) ||
	  (position == 0),
	  ExcInvalidPosition());
				   // first insert the new time step
				   // into the doubly linked list
				   // of timesteps
  if (position == 0)
    {
				       // at the end
      new_timestep->set_next_timestep (0);
      if (timesteps.size() > 0)
	{
	  timesteps.back()->set_next_timestep (new_timestep);
	  new_timestep->set_previous_timestep (timesteps.back());
	}
      else
	new_timestep->set_previous_timestep (0);
    }
  else
    if (position == timesteps[0])
      {
					 // at the beginning
	new_timestep->set_previous_timestep (0);
	if (timesteps.size() > 0)
	  {
	    timesteps[0]->set_previous_timestep (new_timestep);
	    new_timestep->set_next_timestep (timesteps[0]);
	  }
	else
	  new_timestep->set_next_timestep (0);
      }
    else
      {
					 // inner time step
	std::vector<SmartPointer<TimeStepBase,TimeDependent> >::iterator insert_position
	  = std::find(timesteps.begin(), timesteps.end(), position);

	(*(insert_position-1))->set_next_timestep (new_timestep);
	new_timestep->set_previous_timestep (*(insert_position-1));
	new_timestep->set_next_timestep (*insert_position);
	(*insert_position)->set_previous_timestep (new_timestep);
      };

				   // finally enter it into the
				   // array
  timesteps.insert ((position == 0 ?
		     timesteps.end() :
		     std::find(timesteps.begin(), timesteps.end(), position)),
		    new_timestep);
}


void
TimeDependent::add_timestep (TimeStepBase *new_timestep)
{
  insert_timestep (0, new_timestep);
}


void TimeDependent::delete_timestep (const unsigned int position)
{
  Assert (position<timesteps.size(),
	  ExcInvalidPosition());

				   // Remember time step object for
				   // later deletion and unlock
				   // SmartPointer
  TimeStepBase* t = timesteps[position];
  timesteps[position] = 0;
				   // Now delete unsubscribed object
  delete t;

  timesteps.erase (timesteps.begin() + position);

				   // reset "next" pointer of previous
				   // time step if possible
				   //
				   // note that if now position==size,
				   // then we deleted the last time step
  if (position != 0)
    timesteps[position-1]->set_next_timestep ((position<timesteps.size()) ?
					      timesteps[position] :
					      /*null*/SmartPointer<TimeStepBase,TimeDependent>());

				   // same for "previous" pointer of next
				   // time step
  if (position<timesteps.size())
    timesteps[position]->set_previous_timestep ((position!=0) ?
						timesteps[position-1] :
						/*null*/SmartPointer<TimeStepBase,TimeDependent>());
}


void
TimeDependent::solve_primal_problem ()
{
  do_loop (std::mem_fun(&TimeStepBase::init_for_primal_problem),
	   std::mem_fun(&TimeStepBase::solve_primal_problem),
	   timestepping_data_primal,
	   forward);
}


void
TimeDependent::solve_dual_problem ()
{
  do_loop (std::mem_fun(&TimeStepBase::init_for_dual_problem),
	   std::mem_fun(&TimeStepBase::solve_dual_problem),
	   timestepping_data_dual,
	   backward);
}


void
TimeDependent::postprocess ()
{
  do_loop (std::mem_fun(&TimeStepBase::init_for_postprocessing),
	   std::mem_fun(&TimeStepBase::postprocess_timestep),
	   timestepping_data_postprocess,
	   forward);
}



void TimeDependent::start_sweep (const unsigned int s)
{
  sweep_no = s;

				   // reset the number each
				   // time step has, since some time
				   // steps might have been added since
				   // the last time we visited them
				   //
				   // also set the sweep we will
				   // process in the sequel
  for (unsigned int step=0; step<timesteps.size(); ++step)
    {
      timesteps[step]->set_timestep_no (step);
      timesteps[step]->set_sweep_no (sweep_no);
    };

  for (unsigned int step=0; step<timesteps.size(); ++step)
    timesteps[step]->start_sweep ();
}



void TimeDependent::end_sweep (const unsigned int n_threads)
{
  if (DEAL_II_USE_MT && (n_threads > 1))
    {
      const unsigned int stride = timesteps.size() / n_threads;
      Threads::ThreadGroup<> threads;
      void (TimeDependent::*p) (const unsigned int, const unsigned int)
        = &TimeDependent::end_sweep;
      for (unsigned int i=0; i<n_threads; ++i)
        threads += Threads::new_thread (p, *this, i*stride,
					(i == n_threads-1 ?
					 timesteps.size() :
					 (i+1)*stride));
      threads.join_all();
    }
  else
                                     // now do the work
    end_sweep (0, timesteps.size());
}



void TimeDependent::end_sweep (const unsigned int begin,
			       const unsigned int end)
{
  for (unsigned int step=begin; step<end; ++step)
    timesteps[step]->end_sweep ();
}



unsigned int TimeDependent::memory_consumption () const
{
  unsigned int mem = (MemoryConsumption::memory_consumption (timesteps) +
		      MemoryConsumption::memory_consumption (sweep_no) +
		      sizeof(timestepping_data_primal) +
		      sizeof(timestepping_data_dual) +
		      sizeof(timestepping_data_postprocess));
  for (unsigned int i=0; i<timesteps.size(); ++i)
    mem += MemoryConsumption::memory_consumption (*timesteps[i]);

  return mem;
}



/* --------------------------------------------------------------------- */


TimeStepBase::TimeStepBase (const double time) :
		previous_timestep(0),
		next_timestep (0),
		sweep_no (numbers::invalid_unsigned_int),
		timestep_no (numbers::invalid_unsigned_int),
		time (time)
{}



TimeStepBase::~TimeStepBase ()
{}



void
TimeStepBase::wake_up (const unsigned )
{}



void
TimeStepBase::sleep (const unsigned)
{}



void
TimeStepBase::start_sweep ()
{}



void
TimeStepBase::end_sweep ()
{}



void
TimeStepBase::init_for_primal_problem ()
{
  next_action = primal_problem;
}



void
TimeStepBase::init_for_dual_problem ()
{
  next_action = dual_problem;
}



void
TimeStepBase::init_for_postprocessing ()
{
  next_action = postprocess;
}



void
TimeStepBase::solve_dual_problem ()
{
  Assert (false, ExcPureFunctionCalled());
}



void
TimeStepBase::postprocess_timestep ()
{
  Assert (false, ExcPureFunctionCalled());
}



double
TimeStepBase::get_time () const
{
  return time;
}



unsigned int
TimeStepBase::get_timestep_no () const
{
  return timestep_no;
}



double
TimeStepBase::get_backward_timestep () const
{
  Assert (previous_timestep != 0, ExcCantComputeTimestep());
  return time - previous_timestep->time;
}



double
TimeStepBase::get_forward_timestep () const
{
  Assert (next_timestep != 0, ExcCantComputeTimestep());
  return next_timestep->time - time;
}



void
TimeStepBase::set_previous_timestep (const TimeStepBase *previous)
{
  previous_timestep = previous;
}



void
TimeStepBase::set_next_timestep (const TimeStepBase *next)
{
  next_timestep     = next;
}



void
TimeStepBase::set_timestep_no (const unsigned int step_no)
{
  timestep_no = step_no;
}



void
TimeStepBase::set_sweep_no (const unsigned int sweep)
{
  sweep_no = sweep;
}



unsigned int
TimeStepBase::memory_consumption () const
{
				   // only simple data types
  return sizeof(*this);
}


DEAL_II_NAMESPACE_CLOSE
