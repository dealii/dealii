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

// verify that TaskGroup does what we want

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <unistd.h>

#include <base/thread_management.h>


void test (int i) 
{
  deallog << "Task " << i << " starting..." << std::endl;
  sleep (1);
  deallog << "Task " << i << " finished!" << std::endl;
}

  
  

int main()
{
  std::ofstream logfile("task_07/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  Threads::TaskGroup<> tg;
  tg += Threads::new_task (test, 1);
  tg += Threads::new_task (test, 2);

  tg.join_all ();
  
  deallog << "OK" << std::endl;
}
