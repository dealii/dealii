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


// this is a test opposite of thread_validity_05: if we have a function that
// gets an argument *by value* rather than by reference, then make sure that
// it is actually copied when passed to a new thread by reference
//
// deal.II does not guarantee how often an argument is copied, but it has to
// be copied at least once from the calling thread to the stack of the called
// thread

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <deal.II/base/thread_management.h>

struct X
{
  X() : i(0) {}
  X(const X &x) : i(x.i+1) {}
  X &operator= (const X &x)
  {
    i = x.i+1;
    return *this;
  }
  int i;
};


void execute_ref (const X &x)
{
  Assert (x.i == 0, ExcInternalError());
  deallog << __PRETTY_FUNCTION__ << ' ' << x.i << std::endl;
  deallog << "OK" << std::endl;
}

void execute_value (X x)
{
  Assert (x.i > 0, ExcInternalError());
  deallog << __PRETTY_FUNCTION__ << ' ' << (x.i>0 ? "OK" : "not OK")
          << std::endl;
  deallog << "OK" << std::endl;
}


void test ()
{
  {
    X x;
    Threads::Thread<void> t = Threads::spawn (&execute_ref)(x);
    t.join ();
  }
  {
    X x;
    Threads::Thread<void> t = Threads::spawn (&execute_value)(x);
    t.join ();
  }
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog.detach();
  logfile.close();
  unify_pretty_function ("output");
}
