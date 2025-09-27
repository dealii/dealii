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


// verify that thread local storage works as advertised. like _01 but using
// the initialization with an exemplar

#include <deal.II/base/thread_local_storage.h>
#include <deal.II/base/thread_management.h>

#include <atomic>

#include "../tests.h"


struct X
{
  X()
  {
    Assert(false, ExcInternalError());
  };
  X(int n)
  {
    deallog << "Creating" << std::endl;
    Assert(n == 42, ExcInternalError());
  };
  X(const X &)
  {
    deallog << "Copying" << std::endl;
  };
  ~X()
  {
    deallog << "Destroying " << std::endl;
  };
  int i;
};

Threads::ThreadLocalStorage<X> *tls_data;

static std::atomic<int> counter(0);

void
execute(int i)
{
  tls_data->get().i = i;

  // indicate that the TLS object has been
  // accessed
  static Threads::Mutex m;
  {
    std::lock_guard<std::mutex> l(m);
    ++counter;
  }

  // wait in order to make sure that the
  // thread lives longer than the TLS object
  std::this_thread::sleep_for(std::chrono::seconds(5));
}


void
test()
{
  // create a thread local storage object
  X exemplar(42);
  tls_data = new Threads::ThreadLocalStorage<X>(exemplar);

  // create 5 threads and wait for their
  // return. the OS will create 5 individual
  // thread ids, which means that we will
  // create 5 individual thread specific
  // storage locations
  std::vector<std::thread> tg;
  for (unsigned int i = 10; i < 15; ++i)
    tg.emplace_back(execute, i);

  // spin lock until all threads have created
  // their objects
  while (counter != 5)
    ;

  // delete the TLS object. this should also
  // destroy all the objects created so far,
  // namely 5 of them, plus the one copy the
  // TLS object stores; while that happens,
  // we should get output into logstream
  delete tls_data;

  // write something into the output
  // file. this also makes sure that the
  // output file records the order in which
  // this output is generated and in which
  // the TLS objects were destroyed. we want
  // that these objects are destroyed when
  // the TLS object is destroyed, not when
  // the threads exit
  deallog << "Done." << std::endl;

  // now make sure the threads all finish
  for (auto &thread : tg)
    thread.join();

  // at this point, the seventh object will
  // be destroyed, which is the exemplar
  // object local to this function
}



int
main()
{
  initlog();

  test();
}
