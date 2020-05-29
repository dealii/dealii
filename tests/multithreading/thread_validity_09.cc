// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// see if we can detach from threads. before r18272 we used to have a
// bug where detached threads would write their return value into
// released memory. this test releases memory, allocates it again,
// and makes sure nobody writes into it at undue times

#include <deal.II/base/thread_management.h>

#include <atomic>

#include "../tests.h"


Threads::Mutex          mutex;
static std::atomic<int> spin_lock(0);


int
worker()
{
  // wait till the main thread has actually done its work -- we will
  // hang in the 'acquire' line until the main thread releases the
  // mutex. release the mutex again at the end of this function since
  // mutices can only be released on the same thread as they are
  // acquired on.
  mutex.lock();
  deallog << "OK." << std::endl;
  spin_lock = 1;
  mutex.unlock();
  return 42;
}



int
main()
{
  initlog();

  mutex.lock();
  {
    Threads::new_thread(worker);
  }
  std::this_thread::sleep_for(std::chrono::seconds(1));

  const unsigned int sz = 1000000;
  char *             p  = new char[sz];
  for (unsigned int i = 0; i < sz; ++i)
    p[i] = 0;

  // make sure the worker thread can actually start
  mutex.unlock();

  // wait for the worker thread to do its work
  while (spin_lock == 0)
    ;

  for (unsigned int i = 0; i < sz; ++i)
    AssertThrow(p[i] == 0, ExcInternalError());

  delete[] p;
}
