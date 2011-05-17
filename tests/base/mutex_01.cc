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

// verify that mutexes work correctly in MT context

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <unistd.h>

#include <deal.II/base/thread_management.h>


Threads::ThreadMutex mutex;


void test () 
{
                                   // get the mutex, but note that it is first
                                   // held by main() which will therefore
                                   // usually get to write to deallog first
  mutex.acquire();
  deallog << "2" << std::endl;
}

  
int main()
{
  std::ofstream logfile("mutex_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  if (DEAL_II_USE_MT != 0)
    {
      mutex.acquire ();
  
      Threads::Thread<> t = Threads::new_thread (&test);

      sleep (2);
      deallog << "1" << std::endl;
  
      mutex.release ();
      t.join ();

				       // the other thread should now have
				       // acquired the mutex, so release
				       // it again (destroying an acquired
				       // lock invokes undefined behavior
				       // in pthread_mutex_destroy)
      mutex.release ();
    }
  else
    {
				       // make sure the test also
				       // works in non-MT mode
      deallog << "1" << std::endl;
      deallog << "2" << std::endl;
    }
}
