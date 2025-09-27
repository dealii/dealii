// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// the same as timer.cc but also test the new functions last_wall_time(),
// cpu_time() and last_cpu_time().

#include <deal.II/base/timer.h>

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


int
main()
{
  initlog();

  Timer t;
  burn(50);

  const double old_wall_time = t.wall_time();
  AssertThrow(old_wall_time > 0., ExcInternalError());
  const double old_cpu_time = t.cpu_time();
  AssertThrow(old_cpu_time > 0., ExcInternalError());

  burn(50);
  AssertThrow(t.stop() > 0., ExcInternalError());

  AssertThrow(t.wall_time() > old_wall_time, ExcInternalError());
  AssertThrow(t.cpu_time() > old_cpu_time, ExcInternalError());
  AssertThrow(t.last_wall_time() > 0., ExcInternalError());
  AssertThrow(t.last_cpu_time() > 0, ExcInternalError());

  t.reset();
  AssertThrow(t.wall_time() == 0., ExcInternalError());
  AssertThrow(t.cpu_time() == 0., ExcInternalError());

  deallog << "OK" << std::endl;
}
