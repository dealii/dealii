//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// see that we can query a thread object that has never been
// assigned

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <unistd.h>

#include <deal.II/base/thread_management.h>


Threads::Mutex mutex;
int spin_lock = 0;


int worker ()
{
  sleep (1);
  return 42;
}



int main()
{
  std::ofstream logfile("thread_validity_10/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Threads::Thread<int> t;
				   // join non-existing thread
  deallog << (t.valid() ? "true" : "false")
	  << std::endl;

				   // now assign a thread object and
				   // wait for it
  t = Threads::new_thread (worker);
  deallog << (t.valid() ? "true" : "false")
	  << std::endl;
  deallog << "return value = " << t.return_value()
	  << std::endl;
}
