//----------------------------  stack_trace.cc  ---------------------------
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
//----------------------------  stack_trace.cc  ---------------------------


#include "../tests.h"
#include <base/exceptions.h>
#include <base/logstream.h>
#include <base/function_lib.h>
#include <fstream>

// test the stack trace generation when hitting an Assert(...)


void test1a ()
{
  int i = 0;
  Assert (i!=0, ExcMessage ("Gotcha!"));
}


void test1b ()
{
  int i = 0;
  AssertThrow (i!=0, ExcMessage ("Gotcha!"));
}


void test2 ()
{

				   // access an invalid
				   // component. should trigger an
				   // Assert()
  ConstantFunction<2>(1.).value (Point<2>(), 2);
}



int main ()
{
  std::ofstream logfile("stack_trace.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try { test1a (); }
  catch (std::exception &exc) {
    deallog << "  caught exception:" << std::endl
	    << exc.what() << std::endl;
  }

  try { test1b (); }
  catch (std::exception &exc) {
    deallog << "  caught exception:" << std::endl
	    << exc.what() << std::endl;
  }

  try { test2 (); }
  catch (std::exception &exc) {
    deallog << "  caught exception:" << std::endl
	    << exc.what() << std::endl;
  }
}

