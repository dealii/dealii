// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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


// TimerOutput is calling MPI functions internally when deal.II is
// configured with MPI. This creates an error of the form:
// *** The MPI_Allreduce() function was called before MPI_INIT was invoked.
// *** This is disallowed by the MPI standard.
// *** Your MPI job will now abort.

#include "../tests.h"
#include <deal.II/base/timer.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <cmath>
#include <iomanip>

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-1);

  {
    
    // use std::cout so that no output is saved to the logfile, because it
    // is difficult to test (timing)
    TimerOutput t(std::cout, TimerOutput::summary, TimerOutput::cpu_times);

    t.enter_subsection("hi");
    t.leave_subsection("hi");
  }

  deallog << "ok" << std::endl;
}

