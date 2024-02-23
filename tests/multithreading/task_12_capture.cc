// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// verify that we can create tasks both via lambdas and via std::bind
// expressions. this obviously requires C++11
//
// like the _12 test, but using lambda captures


#include <deal.II/base/thread_management.h>

#include "../tests.h"


// return a double, to make sure we correctly identify the return type
// of the expressions used in new_task(...)
double
test(int i)
{
  deallog << "Task " << i << " starting..." << std::endl;
  std::this_thread::sleep_for(std::chrono::seconds(1));
  deallog << "Task " << i << " finished!" << std::endl;

  return 3.141;
}



int
main()
{
  initlog();

  int arg1 = 1;
  int arg2 = 2;

  Threads::TaskGroup<double> tg;
  tg += Threads::new_task([arg1]() // capture arg1 by value
                          { return test(arg1); });
  tg += Threads::new_task([&arg2]() // capture arg2 by reference
                          { return test(arg2); });

  tg.join_all();

  deallog << "OK" << std::endl;

  std::ofstream *out_stream =
    dynamic_cast<std::ofstream *>(&deallog.get_file_stream());
  Assert(out_stream != nullptr, ExcInternalError());
  deallog.detach();
  out_stream->close();
  sort_file_contents("output");
}
