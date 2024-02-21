// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// verify that a change to TimeDependent::end_sweep works as expected:
// TimeStep::end_sweep must be called once for every time step object


#include <deal.II/numerics/time_dependent.h>

#include <algorithm>

#include "../tests.h"


std::vector<bool> end_sweep_flags;


class TimeStep : public TimeStepBase
{
public:
  TimeStep(const unsigned int time_step_number)
    : TimeStepBase(0)
    , time_step_number(time_step_number)
  {}

  virtual void
  end_sweep()
  {
    static Threads::Mutex       mutex;
    std::lock_guard<std::mutex> lock(mutex);
    end_sweep_flags[time_step_number] = true;
  }

  virtual void
  solve_primal_problem()
  {}

private:
  const unsigned int time_step_number;
};


void
test()
{
  // create time steps, more than there are likely threads on current machines
  TimeDependent      td(TimeDependent::TimeSteppingData(0, 0),
                   TimeDependent::TimeSteppingData(0, 0),
                   TimeDependent::TimeSteppingData(0, 0));
  const unsigned int n_time_steps = 10000;
  for (unsigned int i = 0; i < n_time_steps; ++i)
    td.add_timestep(new TimeStep(i));

  end_sweep_flags = std::vector<bool>(n_time_steps, false);
  td.end_sweep();

  // make sure we have called TimeStep::end_sweep once for every time step
  // object
  AssertThrow(end_sweep_flags == std::vector<bool>(n_time_steps, true),
              ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test();

  return 0;
}
