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

// Make sure we can call Threads::Thread::join on objects that haven't even
// been assigned a thread

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <deal.II/base/thread_management.h>

void execute ()
{}


void test () 
{
				   // use a default constructed object
  Threads::Thread<> t;
  deallog << "Before first join()" << std::endl;
  t.join ();
  deallog << "Between join()s" << std::endl;
  t.join ();
  deallog << "After second join()" << std::endl;  
}

  
  

int main()
{
  std::ofstream logfile("thread_validity_12/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
