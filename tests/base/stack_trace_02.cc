//----------------------------  stack_trace_02.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  stack_trace_02.cc  ---------------------------


#include "../tests.h"
#include <base/exceptions.h>
#include <base/logstream.h>
#include <base/function_lib.h>
#include <fstream>

// test the stack trace generation when hitting an Assert(...)


void test1 ()
{
				   // access invalid component
  ConstantFunction<2>(1.).value (Point<2>(), 42);
}

void test2 ()
{
  test1 ();
}


int main ()
{
  std::ofstream logfile("stack_trace_02.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try { test2 (); }
  catch (std::exception &exc) {
    deallog << "  caught exception:" << std::endl
	    << exc.what() << std::endl;
  }
}

