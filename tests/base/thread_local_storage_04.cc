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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// test ThreadLocalStorage::operator= (const T&)

#include <deal.II/base/thread_local_storage.h>
#include <deal.II/base/thread_management.h>

#include "../tests.h"

int counter = 10;

struct X
{
  Threads::ThreadLocalStorage<int> tls_data;

  X()
    : tls_data(42)
  {}

  int
  f()
  {
    // use TLS::operator=
    tls_data = counter++;
    // access TLS data and have it
    // converted to the right data type
    // without the need to call
    // tls_data.get()
    return tls_data;
  }
};


void
test()
{
  X x;
  {
    Threads::Thread<int> t;
    t = Threads::new_thread(&X::f, x);
    AssertThrow(t.return_value() == 10, ExcInternalError());
  }
  {
    Threads::Thread<int> t;
    t = Threads::new_thread(&X::f, x);
    AssertThrow(t.return_value() == 11, ExcInternalError());
  }

  AssertThrow(counter == 12, ExcInternalError());
}



int
main()
{
  initlog();

  test();
  deallog << "OK" << std::endl;
}
