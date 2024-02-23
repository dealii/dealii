// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Verify that we can create tasks both via lambdas that take
// arguments, and then passing these arguments as follow-up arguments
// to new_task(). In other words, this is the equivalent strategy to
// passing a function pointer as the first argument, and the arguments
// to the function following.
//
// task_13 already checks the case of lambda functions that take no
// arguments.


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

  Threads::TaskGroup<double> tg;

  // use variations of ways we can declare lambdas
  tg += Threads::new_task([](int i) -> double { return test(i); }, 1);
  tg += Threads::new_task([](int i) -> double { return (float)test(i); }, 2);

  tg.join_all();

  deallog << "OK" << std::endl;

  std::ofstream *out_stream =
    dynamic_cast<std::ofstream *>(&deallog.get_file_stream());
  Assert(out_stream != nullptr, ExcInternalError());
  deallog.detach();
  out_stream->close();
  sort_file_contents("output");
}
