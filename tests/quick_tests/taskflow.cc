// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test taskflow

#include <deal.II/base/thread_management.h>

#include <taskflow/algorithm/for_each.hpp>
#include <taskflow/taskflow.hpp>

#include <iostream>

using namespace dealii;

void
test1()
{
  auto &executor = MultithreadInfo::get_taskflow_executor();

  std::atomic<unsigned int> counter;
  counter = 0;

  tf::Taskflow taskflow;

  auto incrementor = [&]() { counter.fetch_add(1); };


  tf::Task A = taskflow.emplace(incrementor);
  tf::Task B = taskflow.emplace(incrementor);
  tf::Task C = taskflow.emplace([&]() {
    std::cout << "counter = " << counter << std::endl;
    counter = 0;
  });

  A.precede(C);
  B.precede(C);

  auto p =
    taskflow.for_each_index(1, 11, 1, [&](int idx) { counter.fetch_add(idx); });

  C.precede(p);

  executor.run(taskflow).wait();

  if (counter != 55)
    {
      std::cout << "error: counter is " << counter << " and not 55."
                << std::endl;
      exit(1);
    }

  taskflow.dump(std::cout);
}



int
main()
{
  MultithreadInfo::set_thread_limit();

  std::cout << "taskflow will use "
            << MultithreadInfo::get_taskflow_executor().num_workers()
            << " out of " << MultithreadInfo::n_cores() << " cores."
            << std::endl
            << "MultithreadInfo::n_thread()= " << MultithreadInfo::n_threads()
            << std::endl;

  test1();
}
