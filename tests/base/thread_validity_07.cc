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

// test that if multiple threads are waiting for one single thread, all of the
// waiting ones will be woken up.

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <unistd.h>

#include <base/thread_management.h>


void worker ()
{
  deallog << "Worker thread is starting." << std::endl;
  sleep (3);
  deallog << "Worker thread is finished." << std::endl;
}

Threads::Thread<> worker_thread;

void waiter (int i) 
{
  worker_thread.join ();
  
  deallog << "Waiting thread " << i << " was woken up."
	  << std::endl;
}

  
  

int main()
{
  std::ofstream logfile("thread_validity_07/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  worker_thread = Threads::new_thread (worker);

  Threads::ThreadGroup<> waiter_threads;
  for (unsigned int i=0; i<20; ++i)
    waiter_threads += Threads::new_thread (waiter, i);

  waiter_threads.join_all ();
  deallog << "All waiting threads finished." << std::endl;
}
