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
// this is a variant that makes sure that if there are const and non-const
// member functions that the correct one is called

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <base/thread_management.h>

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
    X(const X&);
    X & operator= (const X&);
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
    Y(const Y&);
    Y & operator= (const Y&);
};


void test2 () 
{
  Y y(42);
  Threads::Thread<void> t = Threads::spawn (y, &Y::execute)();
  t.join ();
}

  
  

int main()
{
  std::ofstream logfile("thread_validity_05/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test1 ();
  test2 ();
}
