// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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


// Check that the time step can be increased when using an embedded method.
// Bug reported by Vaibhav Palkar on the mailing list.

#include <deal.II/base/time_stepping.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

Vector<double>
f(double const t, Vector<double> const &y)
{
  Vector<double> values(y);
  for (unsigned int i = 0; i < values.size(); ++i)
    values[i] = 1.0;

  return values;
}

double
my(double const t)
{
  return t;
}

void
test(TimeStepping::EmbeddedExplicitRungeKutta<Vector<double>> &          solver,
     std::function<Vector<double>(double const, Vector<double> const &)> f,
     std::function<double(double const)>                                 my)
{
  double         initial_time = 0.0, final_time = 1.0;
  double         time_step    = 0.1;
  unsigned int   size         = 1;
  unsigned int   n_time_steps = 0;
  double         time         = initial_time;
  Vector<double> solution(size);
  for (unsigned int i = 0; i < size; ++i)
    solution[i] = my(initial_time);


  while (time < final_time)
    {
      if (time + time_step > final_time)
        time_step = final_time - time;
      time      = solver.evolve_one_time_step(f, time, time_step, solution);
      time_step = solver.get_status().delta_t_guess;
      ++n_time_steps;
    }

  deallog << n_time_steps << std::endl;
}

int
main()
{
  initlog();

  TimeStepping::EmbeddedExplicitRungeKutta<Vector<double>> he(
    TimeStepping::HEUN_EULER);
  test(he, f, my);

  return 0;
}
