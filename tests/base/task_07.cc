// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2017 by the deal.II authors
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


// verify that TaskGroup does what we want

#include <deal.II/base/thread_management.h>

#include "../tests.h"


void
test(int i)
{
  deallog << "Task " << i << " starting..." << std::endl;
  std::this_thread::sleep_for(std::chrono::seconds(1));
  deallog << "Task " << i << " finished!" << std::endl;
}



int
main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  Threads::TaskGroup<> tg;
  tg += Threads::new_task(test, 1);
  tg += Threads::new_task(test, 2);

  tg.join_all();

  deallog << "OK" << std::endl;

  deallog.detach();
  logfile.close();
  sort_file_contents("output");
}
