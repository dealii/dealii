//----------------------------  timer.cc  ---------------------------
//    $Id: timer.cc 23710 2011-05-17 04:50:10Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2010, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  timer.cc  ---------------------------

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
  std::ofstream logfile("timer_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  TimerOutput t(logfile, TimerOutput::summary, TimerOutput::cpu_times);

  t.enter_subsection("hi");
  t.leave_subsection("hi");

  t.print_summary();
  deallog << "ok" << std::endl;
}

