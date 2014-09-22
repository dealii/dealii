// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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


// test that objects that can't be copied aren't copied when passed to a new
// thread by reference
//
// this is a variant that makes sure that if there are const and non-const
// member functions that the correct one is called

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <deal.II/base/thread_management.h>

struct X
{
  X(int i) : i(i) {}
  int i;

  void execute ()
  {
    Assert (false, ExcInternalError());
  }

  void execute () const
  {
    Assert (i == 42, ExcInternalError());
    deallog << "OK" << std::endl;
  }

private:
  X(const X &);
  X &operator= (const X &);
};


void test1 ()
{
  const X x(42);
  Threads::Thread<void> t = Threads::spawn (x, &X::execute)();
  t.join ();
}


struct Y
{
  Y(int i) : i(i) {}
  int i;

  void execute ()
  {
    Assert (i == 42, ExcInternalError());
    deallog << "OK" << std::endl;
  }

  void execute () const
  {
    Assert (false, ExcInternalError());
  }

private:
  Y(const Y &);
  Y &operator= (const Y &);
};


void test2 ()
{
  Y y(42);
  Threads::Thread<void> t = Threads::spawn (y, &Y::execute)();
  t.join ();
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test1 ();
  test2 ();
}
