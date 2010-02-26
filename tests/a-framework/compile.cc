//----------------------------  timer.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2010by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  timer.cc  ---------------------------


// test the testsuite framework. this test is supposed to compile
// successfully but not link

#include "../tests.h"
#include <base/logstream.h>
#include <fstream>
#include <cstdlib>

void function_that_does_not_exist ();


int main ()
{
  std::ofstream logfile("compile/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  function_that_does_not_exist ();
  
  deallog << "OK" << std::endl;
}

