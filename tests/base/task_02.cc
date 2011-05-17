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

// like task_01, but with return value

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <unistd.h>

#include <deal.II/base/thread_management.h>


int test () 
{
  sleep (3);
  return 42;
}

  
  

int main()
{
  std::ofstream logfile("task_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Threads::Task<int> t = Threads::new_task (test);
  Assert (t.return_value() == 42, ExcInternalError());

  deallog << "OK" << std::endl;
}
