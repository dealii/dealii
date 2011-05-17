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

// make sure we can start tasks from individual threads. this requires that a
// task scheduler object is running on each thread we create

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <unistd.h>

#include <deal.II/base/thread_management.h>


void test (int i) 
{
  deallog << "Task " << i << " starting..." << std::endl;
  sleep (1);
  if (i<10)
    {
      Threads::new_task (test, 10+i).join ();
    }
  deallog << "Task " << i << " finished!" << std::endl;
}

  
  

int main()
{
  std::ofstream logfile("task_06/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  Threads::Thread<> t1 = Threads::new_thread (test, 1);
  Threads::Thread<> t2 = Threads::new_thread (test, 2);

  t1.join ();
  t2.join ();
  
  deallog << "OK" << std::endl;
}
