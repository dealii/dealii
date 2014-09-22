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


// It looks like we can't call std::thread::join() twice -- the second call
// produces a std::system_error exception, and this can in fact be verified
// because std::thread::joinable() returns false after the first call to
// join()

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <deal.II/base/thread_management.h>

void execute ()
{}


void test ()
{
  Threads::Thread<> t = Threads::spawn (&execute)();
  deallog << "Before first join()" << std::endl;
  t.join ();
  deallog << "Between join()s" << std::endl;
  t.join ();
  deallog << "After second join()" << std::endl;
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
