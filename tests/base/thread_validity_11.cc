//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2008, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

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
  std::ofstream logfile("thread_validity_11/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
