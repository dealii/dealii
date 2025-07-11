// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/discrete_time.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>

DEAL_II_NAMESPACE_OPEN

namespace
{
  // Helper function that computes the next discrete time, adjusting it if:
  //  - The next time exceeds the end time.
  //  - The next time is smaller but very close to the end time.
  double
  calculate_next_time(const double current_time,
                      const double step_size,
                      const double end_time)
  {
    Assert(step_size >= 0., ExcMessage("Time step size must be non-negative"));
    Assert(end_time >= current_time, ExcInternalError());
    double           next_time          = current_time + step_size;
    constexpr double relative_tolerance = 0.05;
    const double     time_tolerance     = relative_tolerance * step_size;
    if (next_time > end_time - time_tolerance)
      next_time = end_time;
    return next_time;
  }
} // namespace



DiscreteTime::DiscreteTime(const double start_time,
                           const double end_time,
                           const double desired_start_step_size)
  : start_time{start_time}
  , end_time{end_time}
  , current_time{start_time}
  , next_time{calculate_next_time(start_time,
                                  desired_start_step_size,
                                  end_time)}
  , previous_time{start_time}
  , start_step_size{next_time - start_time}
  , step_number{0}
{}



void
DiscreteTime::set_desired_next_step_size(const double next_step_size)
{
  next_time = calculate_next_time(current_time, next_step_size, end_time);
}



void
DiscreteTime::set_next_step_size(const double next_step_size)
{
  Assert(next_step_size > 0,
         ExcMessage("Only positive time step size is allowed."));
  next_time = current_time + next_step_size;
  Assert(
    next_time <= end_time,
    ExcMessage(
      "Time step size is too large. The next time cannot exceed the end time."));
}



void
DiscreteTime::advance_time()
{
  Assert(next_time > current_time,
         ExcMessage("You can't advance time further. "
                    "Either dt==0 or you are at the "
                    "end of the simulation time."));
  const double step_size = get_next_step_size();
  previous_time          = current_time;
  current_time           = next_time;
  ++step_number;
  next_time = calculate_next_time(current_time, step_size, end_time);
}



void
DiscreteTime::restart()
{
  previous_time = start_time;
  current_time  = start_time;
  next_time     = calculate_next_time(current_time, start_step_size, end_time);
  step_number   = 0;
}



std::size_t
DiscreteTime::memory_consumption() const
{
  return (MemoryConsumption::memory_consumption(start_time) +
          MemoryConsumption::memory_consumption(end_time) +
          MemoryConsumption::memory_consumption(current_time) +
          MemoryConsumption::memory_consumption(next_time) +
          MemoryConsumption::memory_consumption(previous_time) +
          MemoryConsumption::memory_consumption(start_step_size) +
          MemoryConsumption::memory_consumption(step_number));
}

DEAL_II_NAMESPACE_CLOSE
