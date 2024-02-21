// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test if ThreadLocalStorage::operator= (ThreadLocalStorage<T> &&t) moves the
// shared ptr to the exemplar

#include <deal.II/base/thread_local_storage.h>
#include <deal.II/base/thread_management.h>

#include "../tests.h"


struct X
{
  X()
  {
    deallog << "Creating" << std::endl;
  };
  X(const X &)
  {
    deallog << "Copying" << std::endl;
  };
  int i;
};

Threads::ThreadLocalStorage<X> tls_data;

void
execute(int i)
{
  tls_data.get().i = i;
}


void
test()
{
  // create a thread local storage object
  X                              exemplar;
  Threads::ThreadLocalStorage<X> temp(exemplar);

  // move assign thread local storage object
  tls_data = std::move(temp);

  // create 5 threads and wait for their
  // return. the OS will create 5 individual
  // thread ids, which means that we will
  // create 5 individual thread specific
  // storage locations
  std::vector<std::thread> tg;
  for (unsigned int i = 10; i < 15; ++i)
    tg.emplace_back(execute, i);

  // now make sure the threads all finish
  for (auto &thread : tg)
    thread.join();
}



int
main()
{
  initlog();

  test();
}
