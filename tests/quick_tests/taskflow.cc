// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

// test taskflow

#include <deal.II/base/thread_management.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <taskflow/taskflow.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

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
    taskflow.parallel_for(1, 11, 1, [&](int idx) { counter.fetch_add(idx); });

  C.precede(p.first);

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
