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


// test that if multiple threads are waiting for one single thread, all of the
// waiting ones will be woken up.

#include <deal.II/base/thread_management.h>

#include "../tests.h"


void
worker()
{
  deallog << "Worker thread is starting." << std::endl;
  std::this_thread::sleep_for(std::chrono::seconds(3));
  deallog << "Worker thread is finished." << std::endl;
}

Threads::Thread<> worker_thread;

void
waiter(int i)
{
  worker_thread.join();

  deallog << "Waiting thread " << i << " was woken up." << std::endl;
}



int
main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  worker_thread = Threads::new_thread(worker);

  Threads::ThreadGroup<> waiter_threads;
  for (unsigned int i = 0; i < 20; ++i)
    waiter_threads += Threads::new_thread(waiter, i);

  waiter_threads.join_all();
  deallog << "All waiting threads finished." << std::endl;

  deallog.detach();
  logfile.close();
  sort_file_contents("output");
}
