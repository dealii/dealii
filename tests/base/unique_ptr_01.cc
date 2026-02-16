// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



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
