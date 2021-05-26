// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2021 by the deal.II authors
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
  Threads::ThreadGroup<> tg;
  for (unsigned int i = 10; i < 15; ++i)
    tg += Threads::new_thread(execute, i);

  // now make sure the threads all finish
  tg.join_all();
}



int
main()
{
  initlog();

  test();
}
