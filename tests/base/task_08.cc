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

// tasks cannot be abandoned, i.e. if you don't explicitly wait for a
// task to finish then the destructor of Task will wait for the task
// to finish

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <unistd.h>

#include <deal.II/base/thread_management.h>


void test (int i) 
{
  sleep (1);
  deallog << "Task " << i << " finished!" << std::endl;
}

  
  

int main()
{
  std::ofstream logfile("task_08/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  Threads::new_task (test, 1);
  
  deallog << "OK" << std::endl;
}
