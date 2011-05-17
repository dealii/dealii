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

// see if we can detach from threads. before r18272 we used to have a
// bug where detached threads would write their return value into
// release memory. the this test releases memory, allocates it again,
// and makes sure nobody writes into it at undue times

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <unistd.h>

#include <deal.II/base/thread_management.h>


Threads::Mutex mutex;
int spin_lock = 0;


int worker ()
{
  mutex.acquire ();
  deallog << "OK." << std::endl;
  spin_lock = 1;
  return 42;
}

  

int main()
{
  std::ofstream logfile("thread_validity_09/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  mutex.acquire ();
  {
    Threads::new_thread (worker);
  }
  sleep (1);

  const unsigned int sz = 1000000;
  char *p = new char[sz];
  for (unsigned int i=0; i<sz; ++i)
    p[i] = 0;

  mutex.release ();

  while (spin_lock == 0)
    ;
  mutex.release ();

  for (unsigned int i=0; i<sz; ++i)
    Assert (p[i] == 0, ExcInternalError());
}
