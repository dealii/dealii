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
  std::ofstream logfile("output");
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
