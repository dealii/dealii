//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

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
    X(const X&);
    X & operator= (const X&);
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
  std::ofstream logfile("thread_validity_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
