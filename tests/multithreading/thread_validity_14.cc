// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2018 by the deal.II authors
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


// See that we can query a task object's return value for cases
// where these returned objects are not copyable but only
// movable. This is relevant, for example, when calling functions that
// return a std::unique_ptr.

#include <deal.II/base/thread_management.h>

#include "../tests.h"


class X
{
public:
  X(int i)
    : value(i)
  {}

  // default constructor
  X()
    : value(13)
  {}

  // delete the copy constructor
  X(const X &) = delete;

  // move constructor. sets the moved-from value to a recognizable
  // value
  X(X &&x)
  {
    value   = x.value;
    x.value = 0;
  }

  // same idea about copy operators
  X &
  operator=(const X &) = delete;

  X &
  operator=(X &&x)
  {
    value   = x.value;
    x.value = 0;

    return *this;
  }


  int value;
};



X
foo()
{
  return X(42);
}



int
main()
{
  initlog();

  Threads::Task<X> t = Threads::new_task(&foo);

  // wait for the thread to return and query its value
  deallog << t.return_value().value << std::endl;

  // we can't copy the return_value() object directly, but we can move
  // it. do so and check that we still get the correct value. then
  // also check that the value of the original return object has been
  // reset in the move constructor/operator
  X x = std::move(t.return_value());
  deallog << x.value << std::endl;

  deallog << t.return_value().value << std::endl;
}
