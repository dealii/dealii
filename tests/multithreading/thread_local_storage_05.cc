// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test ThreadLocalStorage::clear

#include <deal.II/base/thread_local_storage.h>
#include <deal.II/base/thread_management.h>

#include "../tests.h"


struct X
{
  X()
    : i(1){};

  int i;
};

Threads::ThreadLocalStorage<X> tls_data;


void
execute(Threads::Mutex &m)
{
  // check correct default initialization
  bool exists;
  int  i = tls_data.get(exists).i;
  AssertThrow(i == 1, ExcInternalError());
  AssertThrow(exists == false, ExcInternalError());

  // set value
  tls_data.get(exists).i = 2;

  // try again. should have existed this time around
  i = tls_data.get(exists).i;
  AssertThrow(i == 2, ExcInternalError());
  AssertThrow(exists == true, ExcInternalError());

  // wait for the barrier to clear
  m.lock();
  m.unlock();

  // at this point, the tls object should have been cleared and should
  // be back at its original value
  i = tls_data.get(exists).i;
  AssertThrow(i == 1, ExcInternalError());
  AssertThrow(exists == false, ExcInternalError());
}


void
test()
{
  const unsigned int N = 10;
  Threads::Mutex     m[N];

  // start N threads with mutexes locked
  std::thread tg[N];
  for (unsigned int i = 0; i < N; ++i)
    {
      m[i].lock();
      tg[i] = std::thread(execute, std::ref(m[i]));
    }

  // let threads work through their first part
  std::this_thread::sleep_for(std::chrono::seconds(3));

  // then reset the thread local object and release the mutexes so the
  // threads can actually run to an end
  tls_data.clear();
  for (unsigned int i = 0; i < N; ++i)
    m[i].unlock();

  // now make sure the threads all finish
  for (auto &thread : tg)
    thread.join();

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
