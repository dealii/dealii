// ---------------------------------------------------------------------
// $Id$
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


// see if we can detach from threads. before r18272 we used to have a
// bug where detached threads would write their return value into
// released memory. this test releases memory, allocates it again,
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
