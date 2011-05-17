//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// we used to synchronise child tasks with the one that spawned it by
// using a mutex, but this could lead to deadlocks if there are more
// tasks than processors available, see the emails on the mailing
// lists in late nov 2010.
//
// this test verifies that we no longer deadlock

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <unistd.h>

#include <deal.II/base/thread_management.h>
#include <deal.II/base/multithread_info.h>


void test ()
{
  sleep (1);
}


void outer ()
{
				   // wait for some time to make sure
				   // that really all outer tasks have
				   // been started, then start a bunch
				   // of new tasks and wait for them
				   // to finish. it used to be that
				   // the join() function used a mutex
				   // which could only be acquired
				   // once the spawned task has
				   // finished, but since we already
				   // have so many tasks running at
				   // this point, the children never
				   // get to run and so never release
				   // the mutex that we try to acquire
				   // in join()
  sleep (1);

  Threads::new_task (test).join();
}



int main()
{
  std::ofstream logfile("task_09/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  for (unsigned int i=0; i<multithread_info.n_cpus*2; ++i)
    Threads::new_task (outer);

  deallog << "OK" << std::endl;
}
