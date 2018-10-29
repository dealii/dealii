// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2017 by the deal.II authors
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


// see that we can query a thread object that has never been
// assigned

#include <deal.II/base/thread_management.h>

#include "../tests.h"


Threads::Mutex mutex;
int            spin_lock = 0;


int
worker()
{
  std::this_thread::sleep_for(std::chrono::seconds(1));
  return 42;
}



int
main()
{
  initlog();

  Threads::Thread<int> t;
  // join non-existing thread
  deallog << (t.valid() ? "true" : "false") << std::endl;

  // now assign a thread object and
  // wait for it
  t = Threads::new_thread(worker);
  deallog << (t.valid() ? "true" : "false") << std::endl;
  deallog << "return value = " << t.return_value() << std::endl;
}
