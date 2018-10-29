// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II authors
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



#include <memory>

#include "../tests.h"

// counter for how many objects of type X there are
int counter = 0;

struct X
{
  X()
  {
    ++counter;
  }

  X(const X &)
  {
    ++counter;
  }

  ~X()
  {
    --counter;
  }
};



int
main()
{
  initlog();

  // test with plain new/delete
  {
    AssertThrow(counter == 0, ExcInternalError());
    {
      X *p = new X;
      AssertThrow(counter == 1, ExcInternalError());
      delete p;
    }
    AssertThrow(counter == 0, ExcInternalError());
  }

  // test with plain unique_ptr
  {
    AssertThrow(counter == 0, ExcInternalError());
    {
      std::unique_ptr<X> p(new X);
      AssertThrow(counter == 1, ExcInternalError());
    }
    AssertThrow(counter == 0, ExcInternalError());
  }

  // test with plain unique_ptr, but also copy stuff. this only works
  // with move constructors, so test only in C++11 mode
  {
    AssertThrow(counter == 0, ExcInternalError());
    {
      std::unique_ptr<X> p(new X);
      AssertThrow(counter == 1, ExcInternalError());

      std::unique_ptr<X> q = std::move(p);
      AssertThrow(counter == 1, ExcInternalError());
    }
    AssertThrow(counter == 0, ExcInternalError());
  }

  deallog << "OK" << std::endl;
}
