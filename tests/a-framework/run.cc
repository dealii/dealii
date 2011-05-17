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


// test the testsuite framework. this test is supposed to run successfully

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <fstream>


int main ()
{
  std::ofstream logfile("run/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  deallog << "OK" << std::endl;
}

