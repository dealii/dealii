//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// a simple test for task based programming: just create a task and wait for
// its completion. the task sleeps for a bit to make sure that the waiting
// code works alright

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <unistd.h>

#include <base/thread_management.h>


void test () 
{
  sleep (3);
  deallog << "OK" << std::endl;
}

  
int main()
{
  std::ofstream logfile("task_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Threads::Task<> t = Threads::new_task (test);
  t.join ();
}
