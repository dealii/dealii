// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/base/std_cxx11/unique_ptr.h>
#include <fstream>
#include <iomanip>

// counter for how many objects of type X there are
int counter = 0;

struct X
{
  X ()
  {
    ++counter;
  }

  X (const X &)
  {
    ++counter;
  }

  ~X ()
  {
    --counter;
  }
};



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // test with plain new/delete
  {
    AssertThrow (counter == 0, ExcInternalError());
    {
      X *p = new X;
      AssertThrow (counter == 1, ExcInternalError());
      delete p;
    }
    AssertThrow (counter == 0, ExcInternalError());
  }

  // test with plain unique_ptr
  {
    AssertThrow (counter == 0, ExcInternalError());
    {
      std_cxx11::unique_ptr<X> p (new X);
      AssertThrow (counter == 1, ExcInternalError());
    }
    AssertThrow (counter == 0, ExcInternalError());
  }

  // test with plain unique_ptr, but also copy stuff. this only works
  // with move constructors, so test only in C++11 mode
#ifdef DEAL_II_WITH_CXX11
  {
    AssertThrow (counter == 0, ExcInternalError());
    {
      std_cxx11::unique_ptr<X> p (new X);
      AssertThrow (counter == 1, ExcInternalError());

      std_cxx11::unique_ptr<X> q = std::move(p);
      AssertThrow (counter == 1, ExcInternalError());
    }
    AssertThrow (counter == 0, ExcInternalError());
  }
#endif

  deallog << "OK" << std::endl;
}

