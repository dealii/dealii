// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// TimerOutput is calling MPI functions internally when deal.II is
// configured with MPI. This creates an error of the form:
// *** The MPI_Allreduce() function was called before MPI_INIT was invoked.
// *** This is disallowed by the MPI standard.
// *** Your MPI job will now abort.

#include <deal.II/base/timer.h>

#include "../tests.h"

int
main()
{
  initlog();

  {
    // use std::cout so that no output is saved to the logfile, because it
    // is difficult to test (timing)
    TimerOutput t(std::cout, TimerOutput::summary, TimerOutput::cpu_times);

    t.enter_subsection("hi");
    t.leave_subsection("hi");
  }

  deallog << "ok" << std::endl;
}
