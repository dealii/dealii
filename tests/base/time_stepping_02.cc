// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that the time step can be increased when using an embedded method.
// Bug reported by Vaibhav Palkar on the mailing list.

#include <deal.II/base/time_stepping.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

Vector<double>
f(const double t, const Vector<double> &y)
{
  Vector<double> values(y);
  for (unsigned int i = 0; i < values.size(); ++i)
    values[i] = 1.0;

  return values;
}

double
my(const double t)
{
  return t;
}

void
test(TimeStepping::EmbeddedExplicitRungeKutta<Vector<double>>           &solver,
     std::function<Vector<double>(const double, const Vector<double> &)> f,
     std::function<double(const double)>                                 my)
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
