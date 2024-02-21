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
    auto future = std::async(&X::f, x);
    AssertThrow(future.get() == 10, ExcInternalError());
  }
  {
    auto future = std::async(&X::f, x);
    AssertThrow(future.get() == 11, ExcInternalError());
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
