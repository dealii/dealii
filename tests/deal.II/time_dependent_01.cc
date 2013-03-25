//----------------------------  time_dependent_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  time_dependent_01.cc  ---------------------------

// verify that a change to TimeDependent::end_sweep works as expected:
// TimeStep::end_sweep must be called once for every time step object


#include "../tests.h"
#include <deal.II/numerics/time_dependent.h>

#include <fstream>
#include <algorithm>
#include <cmath>


std::ofstream logfile("time_dependent_01/output");


std::vector<bool> end_sweep_flags;


class TimeStep : public TimeStepBase
{
public:
  TimeStep (const unsigned int time_step_number)
  :
  TimeStepBase(0),
  time_step_number (time_step_number)
    {}

  virtual void end_sweep ()
    {
      static Threads::Mutex mutex;
      Threads::Mutex::ScopedLock lock(mutex);
      end_sweep_flags[time_step_number] = true;
    }

  virtual void solve_primal_problem () {}

private:
  const unsigned int time_step_number;
};


void test ()
{
  // create time steps, more than there are likely threads on current machines
  TimeDependent td (TimeDependent::TimeSteppingData(0,0),
		    TimeDependent::TimeSteppingData(0,0),
		    TimeDependent::TimeSteppingData(0,0));
  const unsigned int n_time_steps = 10000;
  for (unsigned int i=0; i<n_time_steps; ++i)
    td.add_timestep (new TimeStep(i));

  end_sweep_flags = std::vector<bool> (n_time_steps, false);
  td.end_sweep ();

  // make sure we have called TimeStep::end_sweep once for every time step object
  Assert (end_sweep_flags == std::vector<bool> (n_time_steps, true),
	  ExcInternalError());

  deallog << "OK" << std::endl;
}


int main ()
{
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  return 0;
}
