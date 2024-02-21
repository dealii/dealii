// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// start tasks from tasks

#include <deal.II/base/thread_management.h>

#include "../tests.h"


void
test(int i)
{
  deallog << "Task " << i << " starting..." << std::endl;

  if (i < 10)
    {
      Threads::Task<> t1 = Threads::new_task(test, 10 * i + 1);
      Threads::Task<> t2 = Threads::new_task(test, 10 * i + 2);

      t1.join();
      t2.join();
    }

  std::this_thread::sleep_for(std::chrono::seconds(1));
  deallog << "Task " << i << " finished!" << std::endl;
}



int
main()
{
  initlog();

  Threads::Task<> t1 = Threads::new_task(test, 1);
  Threads::Task<> t2 = Threads::new_task(test, 2);

  t1.join();
  t2.join();

  deallog << "OK" << std::endl;

  std::ofstream *out_stream =
    dynamic_cast<std::ofstream *>(&deallog.get_file_stream());
  Assert(out_stream != nullptr, ExcInternalError());
  deallog.detach();
  out_stream->close();
  sort_file_contents("output");
}
