//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// see if we can detach from threads

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <unistd.h>

#include <deal.II/base/thread_management.h>


Threads::Mutex mutex;
int spin_lock = 0;


void worker ()
{
  mutex.acquire ();
  deallog << "OK." << std::endl;
  spin_lock = 1;
}

  

int main()
{
  std::ofstream logfile("thread_validity_08/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  mutex.acquire ();
				   // start and abandon the
				   // thread. because we hold the
				   // lock, the started task can not
				   // proceed to print out "OK" until
				   // we release the lock which we do
				   // only after some time to give the
				   // thread a way to start up
				   //
				   // if detachment from threads isn't
				   // possible, then the worker()
				   // function will hang because it
				   // won't be able to acquire the
				   // mutex
  {
    Threads::new_thread (worker);
  }
  sleep (1);

				   // let abandoned thread continue
  mutex.release ();

				   // wait for thread to finish
  while (spin_lock == 0)
    ;

				   // release the lock that the thread
				   // acquired to avoid pthread errors
  mutex.release ();
}
