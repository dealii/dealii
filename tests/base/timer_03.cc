// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



#include <deal.II/base/timer.h>

#include "../tests.h"

// compute the ratio of two measurements and compare to
// the expected value.

void
compare(double t1, double t2, double ratio)
{
  double r = t2 / t1;
  double d = std::fabs(r - ratio) / ratio;

  // relative error < 25%?
  if (d <= .25)
    {
      deallog << "OK" << std::endl;
    }
  else
    {
      deallog << "Ratio " << r << " should be " << ratio << std::endl;
    }
}

void
match(double v1, double v2)
{
  double eps = 1.0e-6;
  if (std::fabs(v1 - v2) < eps)
    {
      deallog << "OK" << std::endl;
    }
  else
    {
      deallog << "Value " << v1 << " should be " << v2 << std::endl;
    }
}

// burn computer time

double s = 0.;
void
burn(unsigned int n)
{
  for (unsigned int i = 0; i < n; ++i)
    {
      for (unsigned int j = 1; j < 100000; ++j)
        {
          s += 1. / j * i;
        }
    }
}


int
main()
{
  initlog();

  Timer       t1, t2;
  TimerOutput tO(std::cout, TimerOutput::summary, TimerOutput::cpu_times);

  tO.enter_subsection("Section1");
  tO.enter_subsection("Section2");
  burn(50);
  tO.leave_subsection("Section2");
  tO.enter_subsection("Section2");
  burn(50);
  tO.leave_subsection("Section2");
  tO.leave_subsection("Section1");

  std::map<std::string, double> cpu_times =
    tO.get_summary_data(TimerOutput::OutputData::total_cpu_time);
  std::map<std::string, double> wall_times =
    tO.get_summary_data(TimerOutput::OutputData::total_wall_time);
  std::map<std::string, double> calls =
    tO.get_summary_data(TimerOutput::OutputData::n_calls);

  match(calls["Section1"], 1.0);
  match(calls["Section2"], 2.0);

  if (cpu_times["Section1"] * cpu_times["Section2"] > 0.)
    deallog << "OK" << std::endl;
  else
    deallog << "ERROR - total cpu time 0" << std::endl;

  if (wall_times["Section1"] * wall_times["Section2"] > 0.)
    deallog << "OK" << std::endl;
  else
    deallog << "ERROR - total wall time 0" << std::endl;

  compare(cpu_times["Section1"], cpu_times["Section2"], 1.);
  compare(wall_times["Section1"], wall_times["Section2"], 1.);
}
