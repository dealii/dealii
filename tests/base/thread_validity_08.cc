// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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


// see if we can detach from threads

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <unistd.h>

#include <deal.II/base/thread_management.h>


Threads::Mutex mutex;
volatile int spin_lock = 0;


void worker ()
{
  // wait for the mutex to make sure the main
  // thread has already moved on. we can immediately
  // release the mutex again.
  mutex.acquire ();
  mutex.release ();
  deallog << "OK." << std::endl;
  spin_lock = 1;
}



int main()
{
  std::ofstream logfile("output");
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
}
