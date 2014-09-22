// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


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
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

#ifdef DEAL_II_WITH_THREADS

  mutex.acquire ();

  Threads::Thread<> t = Threads::new_thread (&test);

  sleep (2);
  deallog << "1" << std::endl;

  mutex.release ();
  t.join ();

  // the other thread should now have acquired the mutex, so release it
  // again (destroying an acquired lock invokes undefined behavior in
  // pthread_mutex_destroy)
  mutex.release ();

#else

  // make sure the test also works in non-MT mode
  deallog << "1" << std::endl;
  deallog << "2" << std::endl;
#endif
}
