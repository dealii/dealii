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

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <deal.II/base/thread_management.h>

struct X
{
  X(int i) : i(i) {}
  int i;
private:
  X(const X &);
  X &operator= (const X &);
};


void execute (const X &x)
{
  Assert (x.i == 42, ExcInternalError());
  deallog << "OK" << std::endl;
}


void test ()
{
  X x(42);
  Threads::Thread<void> t = Threads::spawn (&execute)(x);
  t.join ();
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
