// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test TimerOutput with standard section names

#include <deal.II/base/timer.h>

#include <algorithm>
#include <sstream>

#include "../tests.h"

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

void
test(TimerOutput::OutputType output_type)
{
  std::stringstream ss;
  {
    TimerOutput t(ss, TimerOutput::summary, output_type);

    t.enter_subsection("Hello? Hello? Hello?");
    burn(50);
    t.leave_subsection("Hello? Hello? Hello?");

    t.enter_subsection("Is there anybody in there?");
    burn(50);
    t.leave_subsection("Is there anybody in there?");
  }

  std::string s = ss.str();
  std::replace_if(s.begin(), s.end(), ::isdigit, ' ');
  std::replace_if(
    s.begin(), s.end(), [](char x) { return x == '.'; }, ' ');
  deallog << s << std::endl << std::endl;
}

int
main()
{
  initlog();

  deallog << "cpu_times:" << std::endl;
  test(TimerOutput::cpu_times);
  deallog << "wall_times:" << std::endl;
  test(TimerOutput::wall_times);
  deallog << "cpu_and_wall_times:" << std::endl;
  test(TimerOutput::cpu_and_wall_times);
  deallog << "cpu_and_wall_times_grouped:" << std::endl;
  test(TimerOutput::cpu_and_wall_times_grouped);
}
