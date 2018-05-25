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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// tasks cannot be abandoned, i.e. if you don't explicitly wait for a
// task to finish then the destructor of Task will wait for the task
// to finish

#include <deal.II/base/thread_management.h>

#include <unistd.h>

#include "../tests.h"


void
test(int i)
{
  sleep(1);
  deallog << "Task " << i << " finished!" << std::endl;
}



int
main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  {
    Threads::new_task(test, 1);

    deallog << "OK" << std::endl;
  }

  deallog.detach();
  logfile.close();
  sort_file_contents("output");
}
