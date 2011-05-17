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
//
// this is a variant that makes sure that member functions of objects that
// can't be copied aren't called on copies

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
	Assert (i == 42, ExcInternalError());
	deallog << "OK" << std::endl;
      }
    
  private:
    X(const X&);
    X & operator= (const X&);
};




void test () 
{
  X x(42);
  Threads::Thread<void> t = Threads::spawn (x, &X::execute)();
  t.join ();
}

  
  

int main()
{
  std::ofstream logfile("thread_validity_03/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
