// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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


// verify that a change to TimeDependent::end_sweep works as expected:
// TimeStep::end_sweep must be called once for every time step object


#include "../tests.h"
#include <deal.II/numerics/time_dependent.h>

#include <fstream>
#include <algorithm>
#include <cmath>


std::ofstream logfile("output");


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
