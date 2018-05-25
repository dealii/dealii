// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2017 by the deal.II authors
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


// test ThreadLocalStorage::clear but using a sample object to
// initialize the thread local storage

#include <deal.II/base/thread_local_storage.h>
#include <deal.II/base/thread_management.h>

#include <unistd.h>

#include "../tests.h"


struct X
{
  X() : i(1){};

  int i;
};

X
initializer()
{
  X x;
  x.i = 42;
  return x;
}

X forty_two = initializer();


Threads::ThreadLocalStorage<X> tls_data(forty_two);


void
execute(Threads::Mutex &m)
{
  // check correct default initialization
  bool exists;
  int  i = tls_data.get(exists).i;
  AssertThrow(i == 42, ExcInternalError());
  AssertThrow(exists == false, ExcInternalError());

  // set value
  tls_data.get(exists).i = 2;

  // try again. should have existed this time around
  i = tls_data.get(exists).i;
  AssertThrow(i == 2, ExcInternalError());
  AssertThrow(exists == true, ExcInternalError());

  // wait for the barrier to clear
  m.acquire();
  m.release();

  // at this point, the tls object should have been cleared and should
  // be back at its original value
  i = tls_data.get(exists).i;
  AssertThrow(i == 42, ExcInternalError());
  AssertThrow(exists == false, ExcInternalError());
}


void
test()
{
  const unsigned int N = 10;
  Threads::Mutex     m[N];

  // start N threads with mutices locked
  Threads::ThreadGroup<> tg;
  for (unsigned int i = 0; i < N; ++i)
    {
      m[i].acquire();
      tg += Threads::new_thread(execute, m[i]);
    }

  // let threads work through their first part
  sleep(3);

  // then reset the thread local object and release the mutices so the
  // threads can actually run to an end
  tls_data.clear();
  for (unsigned int i = 0; i < N; ++i)
    m[i].release();

  // now make sure the threads all finish
  tg.join_all();

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
